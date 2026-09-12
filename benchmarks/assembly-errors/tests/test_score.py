"""Analytic adversarial examples independent of the case generator."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'src')]
from benchio import digest,save_json,write_fasta,write_tsv,rc
from score import score
from adapt_agp import read_agp,compose,FIELDS


class ScoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        # Two identical copies have distinct target coordinates. Alignment-based
        # copy reassignment must not hide a missing copy or an extra output copy.
        self.seq='ACGT'*50+'TTGA'*50+'ACGT'*50
        write_fasta(self.p/'target.fa',{'truth':self.seq})
        write_fasta(self.p/'input.fa',{'original':self.seq})
        write_tsv(self.p/'input_to_target.tsv',[dict(input_id='original',input_start=0,input_end=600,
                  target_id='truth',target_start=0,target_end=600,strand='+')],
                  ['input_id','input_start','input_end','target_id','target_start','target_end','strand'])
        save_json(self.p/'truth.json',dict(case_id='test',assembly_sha256=digest(self.p/'input.fa'),
                  target_sha256=digest(self.p/'target.fa'),mapping='input_to_target.tsv',events=[],
                  preservation_controls=[dict(variant_id='copy2',kind='copy',target_id='truth',start=400,end=600)],
                  target_labels={'truth':{'output_group':'A'}},breakpoint_tolerance_bp=10))

    def tearDown(self):self.tmp.cleanup()

    def run_score(self, parts, sequences=None):
        rows=[];out={}
        for oid,spans in parts.items():
            pos=0;seq=''
            for s,e,strand in spans:
                frag=self.seq[s:e];frag=frag if strand=='+' else rc(frag)
                seq+=frag
                rows.append(dict(output_id=oid,output_start=pos,output_end=pos+e-s,input_id='original',input_start=s,input_end=e,strand=strand))
                pos+=e-s
            out[oid]=seq
        write_fasta(self.p/'out.fa',out if sequences is None else sequences)
        write_tsv(self.p/'ledger.tsv',rows,FIELDS)
        return score(self.p,self.p/'input.fa',self.p/'out.fa',self.p/'ledger.tsv')

    def test_reverse_complement_and_rename_are_free(self):
        result=self.run_score({'renamed':[(0,600,'-')]})
        self.assertEqual(result['target_retention'],1)
        self.assertEqual(result['false_breakpoints'],0)
        self.assertEqual(result['target_relation_discrepancy_count'],0)
        self.assertTrue(result['preservation_controls'][0]['continuous'])

    def test_extra_identical_copy_detected(self):
        r=self.run_score({'all':[(0,600,'+')],'duplicate':[(400,600,'+')]})
        self.assertEqual(r['input_duplicated_bp'],200)
        self.assertEqual(r['target_extra_copy_callable_bp'],200)
        self.assertFalse(r['retain_all_exact'])

    def test_missing_identical_copy_detected(self):
        r=self.run_score({'one':[(0,400,'+')]})
        self.assertEqual(r['target_missing_callable_bp'],200)
        self.assertEqual(r['preservation_controls'][0]['retained_bp'],0)

    def test_false_cuts_penalized_even_when_retention_perfect(self):
        r=self.run_score({'left':[(0,450,'+')],'right':[(450,600,'+')]})
        self.assertEqual(r['target_retention'],1)
        self.assertEqual(r['false_breakpoints'],1)
        self.assertFalse(r['preservation_controls'][0]['continuous'])

    def test_wrong_internal_orientation_detected(self):
        r=self.run_score({'bad':[(0,200,'+'),(200,400,'-'),(400,600,'+')]})
        self.assertEqual(r['target_relation_discrepancy_count'],2)

    def test_fabricated_provenance_rejected(self):
        with self.assertRaisesRegex(ValueError,'disagrees'):
            self.run_score({'out':[(0,600,'+')]},{'out':'A'*600})

    def test_forward_missing_interval_gap_error(self):
        r=self.run_score({'deleted':[(0,200,'+'),(400,600,'+')]})
        self.assertEqual(r['target_missing_callable_bp'],200)
        self.assertEqual(r['target_relation_discrepancies'][0]['gap_error_bp'],-200)

    def test_whole_contig_order_is_free(self):
        r=self.run_score({'second':[(300,600,'+')],'first':[(0,300,'+')]})
        self.assertEqual(r['target_relation_discrepancy_count'],0)

    def test_ragtag_correction_agp_is_inverted(self):
        p=self.p/'test.agp';p.write_text('original\t1\t300\t1\tW\tfragment\t1\t300\t+\n')
        rows=read_agp(p,True)
        self.assertEqual(rows[0]['output_id'],'fragment')
        self.assertEqual(rows[0]['input_id'],'original')
        self.assertEqual(rows[0]['input_end'],300)

    def test_composed_reverse_slice(self):
        r=dict(output_id='final',output_start=0,output_end=50,input_id='middle',input_start=25,input_end=75,strand='-')
        p=dict(output_id='middle',output_start=0,output_end=100,input_id='original',input_start=100,input_end=200,strand='-')
        result=compose([r],[p])[0]
        self.assertEqual((result['input_start'],result['input_end'],result['strand']),(125,175,'+'))


if __name__=='__main__':unittest.main()
