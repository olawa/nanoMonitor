import unittest
from PyQt6.QtWidgets import QApplication
import sys

# Mock Worker
class MockBatchWorker:
    def __init__(self, bam, batch, af, ref=None):
        self.bam = bam
        self.batch = batch
        self.running = False
        
    def start(self): 
        self.running = True
        
    def isRunning(self): 
        return self.running
        
    def stop(self): self.running = False
    def wait(self): pass

class MockMonitor:
    def __init__(self):
        self.variant_queue = {}
        self.scanned_amplicons = set()
        self.batch_variant_worker = None
        self.amplicon_variants = {}
        self.reference_path = "ref.fa"
        self.log_msgs = []
        self.table_refreshed = False
        self.file_queue = ["test.bam"] # Mock file queue

    def log(self, msg):
        self.log_msgs.append(msg)
        
    def refresh_table(self):
        self.table_refreshed = True

    # --- Copy-Paste Logic from nanoMonitor.py (simplified) ---
    def process_variant_queue(self):
        if not self.variant_queue: return
        
        if self.batch_variant_worker and self.batch_variant_worker.isRunning():
            return
            
        batch = self.variant_queue.copy()
        self.variant_queue = {}
        self.scanned_amplicons.update(batch.keys())
        
        # Mock Worker Start
        self.batch_variant_worker = MockBatchWorker("test.bam", batch, 0.015)
        self.batch_variant_worker.start()
        
    def on_batch_variant_result(self, name, variants):
        self.amplicon_variants[name] = variants
        self.refresh_table()

    def on_batch_variant_finished(self, success, msg):
        if self.batch_variant_worker:
            self.batch_variant_worker.running = False
        self.process_variant_queue()


class TestBatchLogic(unittest.TestCase):
    def setUp(self):
        self.monitor = MockMonitor()

    def test_queue_processing(self):
        # 1. Add to queue
        self.monitor.variant_queue = {"AMP1": "chr1:100-200"}
        self.monitor.process_variant_queue()
        
        # Check worker started
        self.assertTrue(self.monitor.batch_variant_worker.isRunning())
        self.assertEqual(len(self.monitor.scanned_amplicons), 1)
        
        # 2. Add more while running
        self.monitor.variant_queue = {"AMP2": "chr2:500-600"}
        self.monitor.process_variant_queue()
        
        # Should still be running first batch (worker is mocked as running)
        self.assertEqual(len(self.monitor.scanned_amplicons), 1)
        
        # 3. Finish first batch
        self.monitor.on_batch_variant_finished(True, "Done")
        
        # Should have picked up second batch
        self.assertEqual(len(self.monitor.scanned_amplicons), 2)
        self.assertIn("AMP2", self.monitor.scanned_amplicons)
        
    def test_results_update(self):
        self.monitor.on_batch_variant_result("AMP1", [{"af": 0.5}])
        self.assertIn("AMP1", self.monitor.amplicon_variants)
        self.assertTrue(self.monitor.table_refreshed)

if __name__ == '__main__':
    unittest.main()
