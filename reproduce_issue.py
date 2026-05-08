
import unittest
from unittest.mock import MagicMock
from collections import Counter
import ns_amplicon

class MockRead:
    def __init__(self, name, seq, is_unmapped, ref_name=None, ref_start=None, ref_end=None):
        self.query_name = name
        self.query_sequence = seq
        self.is_unmapped = is_unmapped
        self.reference_name = ref_name
        self.reference_start = ref_start
        self.reference_end = ref_end
        self.query_alignment_start = 0
        self.query_alignment_end = len(seq)

class TestPrimerDiscovery(unittest.TestCase):
    def test_mapped_read_discovery(self):
        # Create a mapped read
        # Sequence: "ATCG" + "A"*100 + "TGCA" (108bp)
        # Start sig: ATCG... (30bp)
        # End sig: ...TGCA (30bp)
        seq = "ATCG" * 8 + "A" * 50 + "TGCA" * 8
        # len = 32 + 50 + 32 = 114
        
        # Mapped read
        read = MockRead("read1", seq, is_unmapped=False, ref_name="chr1", ref_start=100, ref_end=214)
        # Mock query_alignment_start/end to match full sequence for simplicity
        read.query_alignment_start = 0
        read.query_alignment_end = len(seq)
        
        # Run discovery
        counts, valid = ns_amplicon.discover_primers_from_batch([read], start_length=10, end_length=30)
        
        # Check results
        # Should be sequence based now!
        # Start sig: first 30bp of seq
        # End sig: last 30bp of seq
        expected_start = seq[:30]
        expected_end = seq[-30:]
        expected_pair = (expected_start, expected_end)
        
        print(f"Counts keys: {counts.keys()}")
        
        self.assertIn(expected_pair, counts)
        self.assertEqual(counts[expected_pair], 1)
        
        # Ensure it's NOT using coordinates
        coord_pair = ("chr1", 100, 214)
        self.assertNotIn(coord_pair, counts)

if __name__ == '__main__':
    unittest.main()
