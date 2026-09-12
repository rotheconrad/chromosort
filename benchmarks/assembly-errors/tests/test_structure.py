"""Analytic completeness cases specified without the simulation generator."""
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from score_structure import complete_paths, score_structure
from benchio import digest, save_json, write_fasta, write_tsv


def part(tid, a, b, strand='+', out=0):
    return dict(target_id=tid, start=a, end=b, strand=strand,
                output_start=out, output_end=out+b-a)


class CompleteStructureTests(unittest.TestCase):
    def test_exact_ledger_projection_repair_and_fragmentation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_fasta(root / 'input.fa', {'c': 'AAAAGGGG'})
            write_fasta(root / 'target.fa', {'t': 'AAAACCCC'})
            fields = ['input_id', 'input_start', 'input_end', 'target_id', 'target_start', 'target_end', 'strand']
            mapping = [dict(zip(fields, row)) for row in [('c', 0, 4, 't', 0, 4, '+'), ('c', 4, 8, 't', 4, 8, '-')]]
            write_tsv(root / 'mapping.tsv', mapping, fields)
            save_json(root / 'truth.json', dict(case_id='analytic', assembly_sha256=digest(root / 'input.fa'),
                target_sha256=digest(root / 'target.fa'), mapping='mapping.tsv', events=[],
                preservation_controls=[], target_labels={}, breakpoint_tolerance_bp=0))
            ledger_fields = ['output_id', 'output_start', 'output_end', 'input_id', 'input_start', 'input_end', 'strand']
            for fragmented in [False, True]:
                write_fasta(root / 'output.fa', {'x': 'AAAA', 'y': 'CCCC'} if fragmented else {'x': 'AAAACCCC'})
                values = [('x', 0, 4, 'c', 0, 4, '+'),
                          ('y' if fragmented else 'x', 0 if fragmented else 4, 4 if fragmented else 8, 'c', 4, 8, '-')]
                write_tsv(root / 'ledger.tsv', [dict(zip(ledger_fields, row)) for row in values], ledger_fields)
                result = score_structure(root, root / 'input.fa', root / 'output.fa', root / 'ledger.tsv')
                self.assertEqual(result['complete_structure'], not fragmented)
                self.assertEqual(result['frozen_score']['target_relation_discrepancy_count'], 0)
                self.assertTrue(result['frozen_score']['retain_all_exact'])

    def test_split_can_remove_bad_adjacency_without_restoring_target(self):
        result = complete_paths({'a': 10}, {'x': 4, 'y': 6},
                                {'x': [part('a', 0, 4)], 'y': [part('a', 4, 10)]})
        self.assertFalse(result['complete_structure'])
        self.assertEqual(result['complete_target_count'], 0)

    def test_reassembled_adjacent_pieces_are_complete(self):
        result = complete_paths({'a': 10}, {'x': 10},
                                {'x': [part('a', 0, 4), part('a', 4, 10, out=4)]})
        self.assertTrue(result['complete_structure'])

    def test_whole_reverse_and_output_record_order_are_free(self):
        result = complete_paths({'a': 10, 'b': 5}, {'y': 5, 'x': 10},
            {'x': [part('a', 4, 10, '-'), part('a', 0, 4, '-', 6)], 'y': [part('b', 0, 5)]})
        self.assertTrue(result['complete_structure'])

    def test_internal_orientation_and_wrong_order_fail(self):
        for path in [[part('a', 0, 4), part('a', 4, 10, '-', 4)],
                     [part('a', 4, 10), part('a', 0, 4, out=6)]]:
            self.assertFalse(complete_paths({'a': 10}, {'x': 10}, {'x': path})['complete_structure'])

    def test_identical_sequence_cannot_replace_a_different_copy_identity(self):
        result = complete_paths({'copy1': 10, 'copy2': 10}, {'x': 10, 'y': 10},
            {'x': [part('copy1', 0, 10)], 'y': [part('copy1', 0, 10)]})
        self.assertFalse(result['complete_structure'])
        self.assertFalse(result['targets']['copy2'])

    def test_extra_or_missing_outputs_fail(self):
        self.assertFalse(complete_paths({'a': 10}, {'x': 10, 'y': 1},
            {'x': [part('a', 0, 10)], 'y': []})['complete_structure'])
        self.assertFalse(complete_paths({'a': 10, 'b': 10}, {'x': 10},
            {'x': [part('a', 0, 10)]})['complete_structure'])

    def test_added_gap_and_truncation_fail(self):
        self.assertFalse(complete_paths({'a': 10}, {'x': 11},
            {'x': [part('a', 0, 4), part('a', 4, 10, out=5)]})['complete_structure'])
        self.assertFalse(complete_paths({'a': 10}, {'x': 9},
            {'x': [part('a', 0, 9)]})['complete_structure'])

    def test_interchromosomal_concatenation_fails(self):
        self.assertFalse(complete_paths({'a': 5, 'b': 5}, {'x': 10},
            {'x': [part('a', 0, 5), part('b', 0, 5, out=5)]})['complete_structure'])


if __name__ == '__main__':
    unittest.main()
