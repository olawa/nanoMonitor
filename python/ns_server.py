# Filename: ns_server.py
# Created: 2025-12-12
# Purpose: Remote server for NanoStream analysis using ZeroMQ

import sys
import os
import zmq
import pickle
import time
import argparse
import glob
from PyQt6.QtCore import QCoreApplication, QObject, pyqtSlot, QTimer, QThread
from ns_workers import AnalysisWorker

class NanoServer(QObject):
    def __init__(self, rep_port=5555, pub_port=5556, secret=None):
        super().__init__()
        self.rep_port = rep_port
        self.pub_port = pub_port
        self.secret = secret
        
        # ZeroMQ Setup
        self.context = zmq.Context()
        self.rep_socket = self.context.socket(zmq.REP)
        self.rep_socket.bind(f"tcp://*:{self.rep_port}")
        
        self.pub_socket = self.context.socket(zmq.PUB)
        self.pub_socket.bind(f"tcp://*:{self.pub_port}")
        
        print(f"[Server] Listening for commands on port {self.rep_port}")
        print(f"[Server] Publishing updates on port {self.pub_port}")

        self.worker = None

        # Timer to poll ZMQ (integrate ZMQ loop with Qt event loop)
        self.timer = QTimer()
        self.timer.timeout.connect(self.check_sockets)
        self.timer.start(10) # Check every 10ms

    def check_sockets(self):
        try:
            # Non-blocking check
            if self.rep_socket.poll(0):
                msg = self.rep_socket.recv()
                self.handle_command(msg)
        except zmq.ZMQError as e:
            print(f"[Server] ZMQ Error: {e}")

    def handle_command(self, msg_bytes):
        try:
            cmd_data = pickle.loads(msg_bytes)
            
            # Auth Check
            if self.secret:
                token = cmd_data.get('auth_token')
                if token != self.secret:
                    print(f"[Server] Auth Failed. Token: {token}")
                    self.rep_socket.send(pickle.dumps("AUTH_ERROR"))
                    return

            command = cmd_data.get('command')
            print(f"[Server] Received command: {command}")

            if command == 'PING':
                self.rep_socket.send(pickle.dumps("PONG"))
            
            elif command == 'LIST_DIR':
                path = cmd_data.get('path', '.')
                recursive = cmd_data.get('recursive', False)
                try:
                    items = self.list_directory(path, recursive=recursive)
                    self.rep_socket.send(pickle.dumps(items))
                except Exception as e:
                    self.rep_socket.send(pickle.dumps({'error': str(e)}))

            elif command == 'START_ANALYSIS':
                params = cmd_data.get('params')
                self.start_analysis(params)
                self.rep_socket.send_string("OK")
                
            elif command == 'STOP_ANALYSIS':
                self.stop_analysis()
                self.rep_socket.send(pickle.dumps("OK"))
                
            elif command == 'RUN_VARIANT':
                params = cmd_data.get('params', {})
                try:
                    variants = self.run_rsnap_variant(params)
                    self.rep_socket.send(pickle.dumps(variants))
                except Exception as e:
                    self.rep_socket.send(pickle.dumps({'error': str(e)}))
                
            else:
                self.rep_socket.send(pickle.dumps("UNKNOWN_COMMAND"))
                
        except Exception as e:
            print(f"[Server] Error handling command: {e}")
            self.rep_socket.send(pickle.dumps(f"ERROR: {str(e)}"))

    def start_analysis(self, params):
        if self.worker is not None:
            self.stop_analysis()
        
        print(f"[Server] Starting analysis on: {params['bam_file']}")
        
        # Instantiate Worker
        # params should contain: bam_file, mode, config, filters, threads
        self.worker = AnalysisWorker(
            bam_file=params['bam_file'],
            mode=params['mode'],
            config=params['config'],
            filters=params['filters'],
            threads=params.get('threads', 4)
        )
        
        # Connect Signals to Publishing Slots
        self.worker.progress.connect(self.on_progress)
        self.worker.results.connect(self.on_results)
        self.worker.partial_results.connect(self.on_partial_results)
        self.worker.finished_file.connect(self.on_finished_file)
        self.worker.error.connect(self.on_error)
        
        # Start Worker
        self.worker.start()

    def stop_analysis(self):
        if self.worker:
            print("[Server] Stopping active worker...")
            self.worker.stop()
            self.worker.wait() # Wait for thread to finish
            self.worker = None
            print("[Server] Worker stopped.")

    def list_directory(self, path, recursive=False):
        if not path: path = "."
        path = os.path.abspath(os.path.expanduser(path))
        
        if not os.path.exists(path):
            return {'error': "Path does not exist"}
        
        if not os.path.isdir(path):
            return {'error': "Not a directory"}
            
        items = []
        if not recursive:
            # Add parent directory
            parent = os.path.dirname(path)
            items.append({'name': '..', 'type': 'dir', 'path': parent, 'size': 0})
            
            with os.scandir(path) as it:
                for entry in it:
                    try:
                        entry_type = 'dir' if entry.is_dir() else 'file'
                        size = entry.stat().st_size if entry_type == 'file' else 0
                        
                        # Filter for relevant files
                        if entry_type == 'file':
                             if not entry.name.lower().endswith(('.bam', '.fastq', '.fastq.gz', '.fq', '.fq.gz', '.bed', '.gtf', '.gff', '.gz')):
                                 continue
                                 
                        items.append({
                            'name': entry.name,
                            'type': entry_type,
                            'path': entry.path,
                            'size': size
                        })
                    except:
                        pass
        else:
            # Recursive search using os.walk
            for root, dirs, files in os.walk(path):
                for name in files:
                    if name.lower().endswith(('.bam', '.fastq', '.fastq.gz', '.fq', '.fq.gz', '.bed', '.gtf', '.gff', '.gz')):
                        full_path = os.path.join(root, name)
                        size = os.path.getsize(full_path)
                        items.append({
                            'name': name,
                            'type': 'file',
                            'path': full_path,
                            'size': size
                        })

        # Sort: Dirs first, then files
        items.sort(key=lambda x: (x['type'] != 'dir', x['name'].lower()))
        return {'current_path': path, 'items': items}

    def run_rsnap_variant(self, params):
        """Runs rsnap variant calling locally on the server."""
        import subprocess
        import uuid
        import os
        
        bam_path = params.get('bam_path')
        region = params.get('region')
        af = params.get('af', 0.015)
        genome_path = params.get('genome_path')
        
        u_id = str(uuid.uuid4())[:8]
        temp_vcf = f"server_var_{u_id}.vcf"
        
        cmd = ["rsnap", "-b", bam_path, "-p", region, "--variant", "--af", str(af), "--output", temp_vcf]
        if genome_path and os.path.exists(genome_path):
            cmd.extend(["-r", genome_path])
            
        print(f"[Server] Running variant calling: {' '.join(cmd)}")
        res = subprocess.run(cmd, capture_output=True, text=True)
        
        variants = []
        if res.returncode == 0 and os.path.exists(temp_vcf):
            with open(temp_vcf, "r") as f:
                for line in f:
                    if line.startswith("#") or not line.strip(): continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 8:
                        chrom = parts[0]
                        try: pos = int(parts[1])
                        except: continue
                        ref, alt = parts[3], parts[4]
                        info_str = parts[7]
                        
                        af_val = 0.0
                        dp = 0
                        info_dict = {}
                        for item in info_str.split(";"):
                            if "=" in item:
                                k, v = item.split("=", 1)
                                info_dict[k] = v
                        if "AF" in info_dict:
                            try: af_val = float(info_dict["AF"])
                            except: pass
                        if "DP" in info_dict:
                            try: dp = int(info_dict["DP"])
                            except: pass
                        
                        variants.append({
                            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "af": af_val, "depth": dp, "vaf": af_val*100
                        })
            try: os.remove(temp_vcf)
            except: pass
            
        return variants

    # --- Signal Handlers (Publish to ZMQ) ---
    
    @pyqtSlot(int)
    def on_progress(self, val):
        self.publish('PROGRESS', val)

    @pyqtSlot(object)
    def on_results(self, data):
        # Data might be large, pickle it
        self.publish('RESULT', data)

    @pyqtSlot(object)
    def on_partial_results(self, data):
        self.publish('PARTIAL_RESULT', data)

    @pyqtSlot(str)
    def on_finished_file(self, filepath):
        self.publish('FINISHED', filepath)

    @pyqtSlot(str)
    def on_error(self, err_msg):
        self.publish('ERROR', err_msg)

    def publish(self, topic, data):
        # Multipart message: [Topic, Pickled Data]
        try:
            self.pub_socket.send_multipart([
                topic.encode('utf-8'),
                pickle.dumps(data)
            ])
        except Exception as e:
            print(f"[Server] Error publishing {topic}: {e}")

def main():
    parser = argparse.ArgumentParser(description="NanoStream Analysis Server")
    parser.add_argument("--rep-port", type=int, default=5555, help="Port for REP socket (commands)")
    parser.add_argument("--pub-port", type=int, default=5556, help="Port for PUB socket (updates)")
    parser.add_argument("--secret", type=str, default=None, help="Authentication Secret")
    
    args = parser.parse_args()
    
    # Create QApp (required for QThread)
    app = QCoreApplication(sys.argv)
    
    server = NanoServer(rep_port=args.rep_port, pub_port=args.pub_port, secret=args.secret)
    
    print("[Server] Ready. Press Ctrl+C to exit.")
    
    try:
        sys.exit(app.exec())
    except KeyboardInterrupt:
        print("[Server] Shutting down...")

if __name__ == "__main__":
    main()
