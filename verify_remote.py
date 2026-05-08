import sys
import os
import subprocess
import time
from PyQt6.QtCore import QCoreApplication, QTimer, pyqtSlot

# Adjust path to import modules
sys.path.append(os.getcwd())
from ns_workers import RemoteAnalysisWorker

SERVER_PORT = 5555
BAM_FILE = "/Users/olwal516/dev/Genomics_Suite/apps_python/nanoStream/bam_test/CGU_2024_05_M22_251003.sup-dup_v5.2.bam.duplex.fastq.gz.hg38.lrhq.yY.bam"

def main():
    app = QCoreApplication(sys.argv)
    
    print("Starting Server...")
    server_process = subprocess.Popen(
        [sys.executable, "ns_server.py", "--rep-port", str(SERVER_PORT), "--pub-port", str(SERVER_PORT + 1), "--secret", "verify_secret"],
        stdout=sys.stdout,
        stderr=sys.stderr
    )
    
    # Allow server to start
    time.sleep(2)
    
    # Use absolute path to dummy file
    BAM_FILE = os.path.abspath("verify_test.bam")

    print("Starting Client Worker...")
    worker = RemoteAnalysisWorker(
        f"tcp://127.0.0.1:{SERVER_PORT}",
        BAM_FILE,
        "Amplicon",
        {"qc_only": True}, # Config
        {"min_qs": 0, "min_len": 0}, # Filters
        threads=4,
        secret="verify_secret"
    )
    
    @pyqtSlot(int)
    def on_progress(val):
        print(f"CLIENT: Progress {val}")

    @pyqtSlot(object)
    def on_results(data):
        print(f"CLIENT: Received Results (Type: {type(data)})")
        # print(data)

    @pyqtSlot(object)
    def on_partial(data):
        print("CLIENT: Partial Result Received")

    @pyqtSlot(str)
    def on_finished(path):
        print(f"CLIENT: Finished {path}")
        app.quit()

    @pyqtSlot(str)
    def on_error(msg):
        print(f"CLIENT: Error: {msg}")
        app.quit()

    worker.progress.connect(on_progress)
    worker.results.connect(on_results)
    worker.partial_results.connect(on_partial)
    worker.finished_file.connect(on_finished)
    worker.error.connect(on_error)
    
    worker.start()
    
    # Timeout
    QTimer.singleShot(20000, app.quit) # 20s timeout
    
    print("Running Event Loop...")
    app.exec()
    
    print("Stopping Server...")
    worker.stop()
    server_process.terminate()
    server_process.wait()
    print("Done.")

if __name__ == "__main__":
    main()
