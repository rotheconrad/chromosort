"""Analytic proposal geometries, independent of sequence generator/tool outputs."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from score_review_proposals import score_review


def event(name,a,b,**extra):
    return dict(event_id=name,contig='c',shape='point' if a==b else 'interval',start=a,end=b,
                kind='inspect',queue=True,decision='pending',**extra)


def truth(errors=(),controls=()):
    return dict(case_id='analytic',contig_lengths={'c':10000},errors=list(errors),biological_controls=list(controls))


def error(name,a,b):
    return dict(error_id=name,contig='c',shape='point' if a==b else 'interval',start=a,end=b)


class ReviewProposalTests(unittest.TestCase):
    def test_plan_counts_once_matches_point_and_can_burden_biology(self):
        t=truth([error('e',5000,5000)],[dict(control_id='v',contig='c',start=7000,end=8000)])
        plan=event('plan',5000,7500);plan.update(shape='plan',positions=[5000,7500])
        result=score_review(t,[plan])
        self.assertEqual(result['queue_size'],1)
        self.assertEqual(result['covered_errors'],1)
        self.assertEqual(result['biological_burden_events'],1)
        self.assertEqual(result['clean_burden_events'],0)
        self.assertEqual(score_review(truth([error('i',5000,7500)]),[plan])['covered_errors'],0)

    def test_orientation_requires_both_edges_on_one_interval(self):
        t=truth([error('e',2000,4000)])
        value=score_review(t,[event('a',2000,2000),event('b',4000,4000)])
        self.assertEqual(value['covered_errors'],0)
        self.assertEqual(value['covered_error_edges'],2)
        self.assertEqual(score_review(t,[event('c',2000,4000)])['covered_errors'],1)

    def test_tolerance_is_inclusive_and_not_overlap_only(self):
        t=truth([error('e',2000,4000)])
        self.assertEqual(score_review(t,[event('a',1500,4500)])['covered_errors'],1)
        self.assertEqual(score_review(t,[event('a',1499,4500)])['covered_errors'],0)
        self.assertEqual(score_review(t,[event('a',0,10000)])['covered_errors'],0)

    def test_one_to_one_matching_and_cross_frame_redundancy(self):
        t=truth([error('e',2000,4000)])
        result=score_review(t,[event('a',2000,4000,reference='one'),event('b',2001,4001,reference='other')])
        self.assertEqual(result['covered_errors'],1)
        self.assertEqual(result['redundant_geometric_events'],1)

    def test_biology_burden_is_not_error_coverage(self):
        t=truth(controls=[dict(control_id='v',contig='c',start=5000,end=7000)])
        result=score_review(t,[event('a',5200,6800)])
        self.assertIsNone(result['error_proposal_coverage'])
        self.assertEqual(result['biological_burden_events'],1)
        self.assertEqual(result['clean_burden_events'],0)

    def test_clean_and_whole_record_control(self):
        t=truth();t['whole_record_orientation_control']=True
        result=score_review(t,[event('a',0,10000)])
        self.assertEqual(result['clean_burden_events'],1)
        self.assertTrue(result['whole_record_orientation_control'])

    def test_acknowledgments_not_actionable_and_duplicates_rejected(self):
        e=event('keep',0,0);e['queue']=False
        self.assertEqual(score_review(truth(),[e])['queue_size'],0)
        with self.assertRaises(ValueError):score_review(truth(),[e,e])

    def test_boundary_event_matches_point_error(self):
        result=score_review(truth([error('e',5000,5000)]),[event('a',5100,5100)])
        self.assertEqual(result['covered_errors'],1)
        self.assertEqual(result['clean_burden_events'],0)

    def test_out_of_bounds_interval_rejected(self):
        with self.assertRaises(ValueError):score_review(truth(),[event('a',9999,10001)])


if __name__=='__main__':unittest.main()
