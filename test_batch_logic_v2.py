import unittest
from PyQt6.QtWidgets import QApplication
import sys
import os

# Mock Workers
class MockWorkerWithFile:
    def __init__(self, path):
        self.bam_file = path

class MockWorkerWithPath:
    def __init__(self, path):
        self.bam_path = path

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
        self.worker_thread = None
        self.amplicon_variants = {}
        self.reference_path = "ref.fa"
        self.log_msgs = []
        self.file_queue = []

    def log(self, msg):
        self.log_msgs.append(msg)
        
    def refresh_table(self):
        pass

    # --- Simulated Logic from nanoMonitor.py ---
    def process_variant_queue(self):
        if not self.variant_queue: return
        
        if self.batch_variant_worker and self.batch_variant_worker.isRunning():
            return
            
        # BAM Resolution Logic
        bam_path = None
        if hasattr(self, 'worker_thread') and self.worker_thread:
            if hasattr(self.worker_thread, 'bam_file'):
                bam_path = self.worker_thread.bam_file
            elif hasattr(self.worker_thread, 'bam_path'):
                bam_path = self.worker_thread.bam_path
        
        if not bam_path and self.file_queue:
             bam_path = self.file_queue[0]
             
        if not bam_path:
            self.log("Skipping variant scan: No valid BAM found.")
            return

        batch = self.variant_queue.copy()
        self.variant_queue = {}
        self.scanned_amplicons.update(batch.keys())
        
        self.batch_variant_worker = MockBatchWorker(bam_path, batch, 0.015)
        self.batch_variant_worker.start()

    def on_results_simulated(self, res):
        # Logic from on_results regarding detection
        new_candidates = {}
        for name, data in res.get("amplicons", {}).items():
            region = data.get("region")
            # Fallback for Nanoparse which provides raw coords
            if not region:
                c = data.get("chrom")
                s = data.get("start")
                e = data.get("end")
                if c and s is not None and e is not None:
                    region = f"{c}:{s}-{e}"
            
            if region and name not in self.scanned_amplicons:
                new_candidates[name] = region
        
        if new_candidates:
            self.variant_queue.update(new_candidates)
            self.process_variant_queue()


class TestBatchLogicUpdated(unittest.TestCase):
    def setUp(self):
        self.monitor = MockMonitor()

    def test_bam_resolution_bam_file(self):
        # Case 1: Worker has bam_file (e.g. NanoparseWorker)
        self.monitor.worker_thread = MockWorkerWithFile("nanoparse.bam")
        self.monitor.variant_queue = {"AMP1": "region1"}
        self.monitor.process_variant_queue()
        
        self.assertIsNotNone(self.monitor.batch_variant_worker)
        self.assertEqual(self.monitor.batch_variant_worker.bam, "nanoparse.bam")

    def test_bam_resolution_bam_path(self):
        # Case 2: Worker has bam_path
        self.monitor.worker_thread = MockWorkerWithPath("analysis.bam")
        self.monitor.variant_queue = {"AMP1": "region1"}
        self.monitor.process_variant_queue()
        
        self.assertIsNotNone(self.monitor.batch_variant_worker)
        self.assertEqual(self.monitor.batch_variant_worker.bam, "analysis.bam")

    def test_nanoparse_region_construction(self):
        # Setup worker first so BAM check passes
        self.monitor.worker_thread = MockWorkerWithFile("test.bam")
        
        # Case 3: Nanoparse result (no region string, just coords)
        res = {
            "amplicons": {
                "AMP_Nano": {
                    "chrom": "chr1", "start": 100, "end": 200, "count": 50
                }
            }
        }
        self.monitor.on_results_simulated(res)
        
        self.assertIn("AMP_Nano", self.monitor.scanned_amplicons)
        self.assertIsNotNone(self.monitor.batch_variant_worker)
        self.assertEqual(self.monitor.batch_variant_worker.batch["AMP_Nano"], "chr1:100-200")

if __name__ == '__main__':
    unittest.main()
