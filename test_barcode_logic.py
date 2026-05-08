import os
import re
import unittest

class MockMonitor:
    def __init__(self):
        self.detected_barcodes = set()
        self.barcode_stats = {}
    
    def extract_barcode(self, filepath):
        """Same logic as nanoMonitor.py"""
        # 1. Check parent folder
        folder = os.path.basename(os.path.dirname(filepath))
        m = re.search(r'(barcode\d+|BC\d+)', folder, re.IGNORECASE)
        if m: return m.group(1).lower()
        
        # 2. Check filename
        filename = os.path.basename(filepath)
        m = re.search(r'(barcode\d+|BC\d+)', filename, re.IGNORECASE)
        if m: return m.group(1).lower()
        
        return "unknown"

    def aggregate_stats(self, bc, res):
        """Logic from on_results"""
        if bc != "unknown":
            if bc not in self.barcode_stats:
                self.barcode_stats[bc] = {}
            for name, data in res.get("amplicons", {}).items():
                if name not in self.barcode_stats[bc]:
                    self.barcode_stats[bc][name] = {"count":0}
                self.barcode_stats[bc][name]["count"] += data["count"]

class TestBarcodeLogic(unittest.TestCase):
    def setUp(self):
        self.monitor = MockMonitor()

    def test_extract_barcode_folder(self):
        self.assertEqual(self.monitor.extract_barcode("/data/barcode01/read.bam"), "barcode01")
        self.assertEqual(self.monitor.extract_barcode("/data/BC02/read.bam"), "bc02")
        self.assertEqual(self.monitor.extract_barcode("barcode99/file.fastq"), "barcode99")

    def test_extract_barcode_filename(self):
        self.assertEqual(self.monitor.extract_barcode("/data/reads/read_barcode05.fastq"), "barcode05")
        self.assertEqual(self.monitor.extract_barcode("BC12_batch1.bam"), "bc12")

    def test_unknown(self):
        self.assertEqual(self.monitor.extract_barcode("/data/reads/file.bam"), "unknown")

    def test_aggregation(self):
        bc = "barcode01"
        res = {"amplicons": {"AMP1": {"count": 10}, "AMP2": {"count": 5}}}
        self.monitor.aggregate_stats(bc, res)
        
        self.assertEqual(self.monitor.barcode_stats[bc]["AMP1"]["count"], 10)
        
        # Add more
        res2 = {"amplicons": {"AMP1": {"count": 2}}}
        self.monitor.aggregate_stats(bc, res2)
        self.assertEqual(self.monitor.barcode_stats[bc]["AMP1"]["count"], 12)

if __name__ == '__main__':
    unittest.main()
