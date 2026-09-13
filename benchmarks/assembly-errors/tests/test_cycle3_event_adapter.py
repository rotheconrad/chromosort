"""Small independent schema-adapter checks; no benchmark/candidate output fixtures."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from normalize_cycle3_events import normalize
from benchio import digest,save_json


class EventAdapterTests(unittest.TestCase):
    def fixture(self,root):
        case=root/'Dtest';scan=root/'scan';case.mkdir();scan.mkdir()
        sha='a'*64
        save_json(case/'manifest.json',dict(assembly=dict(sha256=sha)))
        row=dict(event_id='e',target='c',action='inspect',position=3000,decision='pending',assembly_sha256=sha,reason='aligned review')
        save_json(scan/'scan.json',dict(assembly=dict(sha256=sha),manifest=dict(sha256=digest(case/'manifest.json')),proposals=[row]))
        event=dict(event_id='e',target='c',kind='orientation_discrepancy',start0=2000,end0=4000,positions=[2000,4000],decision='pending',action='inspect',assembly_sha256=sha,
            supporting_block_ids=['b'],opposing_block_ids=[],overlapping_block_ids=[],implicit_edits=False,reasons=['visible'],loss_stages=[],coordinate_evidence=[dict(position=2000),dict(position=4000)])
        save_json(scan/'review_proposals.json',dict(schema='chromosort-review-proposals-v1',assembly_sha256=sha,manifest_sha256=digest(case/'manifest.json'),
            inputs=[],parameters={},blocks=[dict(block_id='b')],contigs=[],events=[event]))
        return case,scan

    def test_structured_event_joins_decision_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            case,scan=self.fixture(Path(temporary));result=normalize(case,scan)
            self.assertEqual(len(result['events']),1)
            self.assertEqual(result['events'][0]['shape'],'interval')
            self.assertEqual(result['events'][0]['positions'],[2000,4000])

    def test_wrong_coordinate_evidence_or_unknown_block_rejected(self):
        for mutation in ['position','block']:
            with tempfile.TemporaryDirectory() as temporary:
                case,scan=self.fixture(Path(temporary));path=scan/'review_proposals.json';data=json.loads(path.read_text())
                if mutation=='position':data['events'][0]['coordinate_evidence'][1]['position']=3999
                else:data['events'][0]['supporting_block_ids']=['unknown']
                save_json(path,data)
                with self.assertRaises(AssertionError):normalize(case,scan)

    def test_alternative_plan_is_one_event_with_explicit_points(self):
        with tempfile.TemporaryDirectory() as temporary:
            case,scan=self.fixture(Path(temporary));path=scan/'review_proposals.json';data=json.loads(path.read_text())
            data['events'][0].update(kind='alternative_segmentation')
            save_json(path,data);result=normalize(case,scan)
            self.assertEqual(len(result['events']),1)
            self.assertEqual(result['events'][0]['shape'],'plan')


if __name__=='__main__':unittest.main()
