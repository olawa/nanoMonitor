# Filename: ns_workers.py
# Created: 2025-11-21 20:10 CET

from PyQt6.QtCore import QThread, pyqtSignal, QObject
import os
import time
import glob
import zmq
import pickle
import ns_amplicon
import ns_rna
import ns_core 

# --- DIRECTORY WATCHER (Unchanged) ---
class DirectoryWatcher(QThread):
    new_files_found = pyqtSignal(list)
    
    def __init__(self, directory, interval=5):
        super().__init__()
        self.directory = directory
        self.extension = "*.bam"
        self.interval = interval
        self.known_files = set()
        self.running = True

    def run(self):
        if self.directory and os.path.exists(self.directory):
            # Watch for BAM and FASTQ
            extensions = ["*.bam", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
            initial_files = []
            for ext in extensions:
                # Recursive search for barcode folders
                initial_files.extend(glob.glob(os.path.join(self.directory, "**", ext), recursive=True))
            
            self.known_files = set(initial_files)
            if self.known_files:
                self.new_files_found.emit(list(self.known_files))
        
        while self.running:
            time.sleep(self.interval)
            if self.directory and os.path.exists(self.directory):
                current_files = set()
                extensions = ["*.bam", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
                for ext in extensions:
                    current_files.update(glob.glob(os.path.join(self.directory, "**", ext), recursive=True))
                
                new_files = list(current_files - self.known_files)
                
                if new_files:
                    self.known_files.update(new_files)
                    self.new_files_found.emit(new_files)
    
    def stop(self):
        self.running = False

# --- ANALYSIS WORKER (Updated with threads argument) ---
class AnalysisWorker(QThread):
    progress = pyqtSignal(int) 
    results = pyqtSignal(object)
    partial_results = pyqtSignal(object)
    finished_file = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, bam_file, mode, config, filters, threads=8, primer_tolerance=0, collect_metadata=True): # Added threads, tolerance
        super().__init__()
        self.bam_file = bam_file
        self.mode = mode
        self.config = config
        self.filters = filters
        self.threads = threads 
        self.primer_tolerance = primer_tolerance
        self.collect_metadata = collect_metadata
        self.is_running = True

    def run(self):
        try:
            res = None
            
            all_meta = []
            def collect_meta(partial):
                if self.collect_metadata and "metadata" in partial:
                    all_meta.extend(partial["metadata"])
                if self.collect_metadata:
                    self.partial_results.emit(partial)

            if self.mode == "Amplicon":
                 res = ns_amplicon.run_analysis(
                    self.bam_file, 
                    self.config.get("primers"), 
                    self.config.get("genes"), 
                    self.filters,
                    self.threads, 
                    lambda p: self.progress.emit(p),
                    collect_meta,
                    self.isInterruptionRequested,
                    qc_only=self.config.get("qc_only", False),
                    primer_tolerance=self.primer_tolerance
                )
            elif self.mode == "RNA-Seq":
                 res = ns_rna.run_analysis(
                    self.bam_file,
                    self.config.get("genes"),
                    self.filters,
                    self.threads, 
                    lambda p: self.progress.emit(p),
                    collect_meta,
                    self.isInterruptionRequested
                )

            # Calculate Pore Stats
            if all_meta:
                try:
                    from ns_core import calculate_pore_stats
                    pore_stats = calculate_pore_stats(all_meta)
                    if pore_stats:
                        if res is None: res = {}
                        res['pore_stats'] = pore_stats
                except Exception as e:
                    print(f"Pore analysis failed: {e}")
            
            if res:
                self.results.emit(res)
            
            self.finished_file.emit(self.bam_file)

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.error.emit(f"AnalysisWorker error for {os.path.basename(self.bam_file)}: {e}")
            
            if res:
                self.results.emit(res)
            
            self.finished_file.emit(self.bam_file)

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.error.emit(f"AnalysisWorker error for {os.path.basename(self.bam_file)}: {e}")
        finally:
            self.is_running = False
            
    def stop(self):
        self.is_running = False

class RemoteAnalysisWorker(QObject): # Note: Inherits QObject, not QThread, but usage is similar
    # Include all signals from AnalysisWorker
    progress = pyqtSignal(int)
    results = pyqtSignal(object)
    partial_results = pyqtSignal(object)
    partial_results = pyqtSignal(object)
    finished_file = pyqtSignal(str)
    finished = pyqtSignal() # Added to match AnalysisWorker interface
    error = pyqtSignal(str)
    
    def __init__(self, server_address, bam_file, mode, config, filters, threads=8, secret=None):
        super().__init__()
        self.server_address = server_address
        self.bam_file = bam_file
        self.mode = mode
        self.config = config
        self.filters = filters
        self.threads = threads
        self.secret = secret
        
        self.context = zmq.Context()
        self.req_socket = None # Control socket
        self.sub_socket = None # Update socket
        self.running = False
        self.listen_thread = None

    def start(self):
        # We need a thread to listen to the SUB socket so we don't block the GUI
        self.running = True
        self.listen_thread = QThread()
        self.worker_thread = RemoteListener(self.server_address, self.bam_file, self.mode, self.config, self.filters, self.threads, self.secret)
        self.worker_thread.moveToThread(self.listen_thread)
        self.worker_thread.moveToThread(self.listen_thread)
        
        # Connect signals
        self.worker_thread.progress.connect(self.progress)
        self.worker_thread.results.connect(self.results)
        self.worker_thread.partial_results.connect(self.partial_results)
        self.worker_thread.finished_file.connect(self.finished_file)
        self.worker_thread.finished_file.connect(lambda _: self.finished.emit()) # Emit finished when file finishes
        self.worker_thread.error.connect(self.error)
        
        self.listen_thread.started.connect(self.worker_thread.run)
        self.listen_thread.start()

    def stop(self):
        if self.worker_thread:
            self.worker_thread.stop()
            self.listen_thread.quit()
            self.listen_thread.wait()

    def terminate(self):
        """Alias for stop() to match QThread interface."""
        self.stop()

class RemoteListener(QObject):
    progress = pyqtSignal(int)
    results = pyqtSignal(object)
    partial_results = pyqtSignal(object)
    finished_file = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, server_address, bam_file, mode, config, filters, threads, secret=None):
        super().__init__()
        self.server_address = server_address # e.g. "tcp://127.0.0.1:5555"
        self.secret = secret
        self.params = {
            'bam_file': bam_file,
            'mode': mode,
            'config': config,
            'filters': filters,
            'threads': threads
        }
        self.running = True

    def run(self):
        context = zmq.Context()
        
        # 1. Setup Request Socket to Send Command
        # Assumes server_address is the REP port (e.g. 5555). 
        # PUB port is assumed to be REP + 1 (5556) or we need to pass it.
        # For simplicity, let's parse the port or assume convention.
        # Let's assume server_address points to REP port.
        
        req_socket = context.socket(zmq.REQ)
        req_socket.connect(self.server_address)
        
        # Derive PUB address
        base = self.server_address.rsplit(':', 1)[0]
        port = int(self.server_address.rsplit(':', 1)[1])
        pub_address = f"{base}:{port + 1}"
        
        # 2. Setup Subscriber Socket
        sub_socket = context.socket(zmq.SUB)
        sub_socket.connect(pub_address)
        sub_socket.setsockopt_string(zmq.SUBSCRIBE, "") # Subscribe to all
        
        # 3. Send Start Command
        cmd = {'command': 'START_ANALYSIS', 'params': self.params}
        if self.secret: cmd['auth_token'] = self.secret
        
        req_socket.send(pickle.dumps(cmd))
        
        # Wait for Ack
        ack = req_socket.recv_string()
        if ack != "OK":
            self.error.emit(f"Server rejected request: {ack}")
            return

        # 4. Listen loop
        poller = zmq.Poller()
        poller.register(sub_socket, zmq.POLLIN)
        
        while self.running:
            try:
                socks = dict(poller.poll(100)) # 100ms timeout
                if sub_socket in socks and socks[sub_socket] == zmq.POLLIN:
                    # Receive multipart: [Topic, Data]
                    msg = sub_socket.recv_multipart()
                    topic = msg[0].decode('utf-8')
                    data = pickle.loads(msg[1])
                    
                    if topic == 'PROGRESS':
                        self.progress.emit(data)
                    elif topic == 'RESULT':
                        self.results.emit(data)
                    elif topic == 'PARTIAL_RESULT':
                        self.partial_results.emit(data)
                    elif topic == 'FINISHED':
                        self.finished_file.emit(data)
                    elif topic == 'ERROR':
                        self.error.emit(data)
                        
            except zmq.ZMQError as e:
                pass # Context terminated or other benign error on shutdown

        # Cleanup will happen when thread exits or object destroyed

    def stop(self):
        self.running = False
        # Send STOP command
        try:
            ctx = zmq.Context()
            sock = ctx.socket(zmq.REQ)
            sock.connect(self.server_address)
            cmd = {'command': 'STOP_ANALYSIS'}
            if self.secret: cmd['auth_token'] = self.secret
            
            sock.send(pickle.dumps(cmd))
            sock.recv_string() # Wait for ACK
            sock.close()
        except:
            pass
        # self.requestInterruption() # Not a QThread

class RemoteClient:
    """Synchronous Client for simple commands like LIST_DIR"""
    def __init__(self, address, secret=None):
        self.address = address
        self.secret = secret
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.REQ)
        self.socket.connect(address)
        
    def list_dir(self, path, recursive=False):
        cmd = {'command': 'LIST_DIR', 'path': path, 'recursive': recursive}
        if self.secret: cmd['auth_token'] = self.secret
        try:
            self.socket.send(pickle.dumps(cmd))
            
            # Use Poller for timeout
            poller = zmq.Poller()
            poller.register(self.socket, zmq.POLLIN)
            if poller.poll(2000): # 2s timeout
                resp = self.socket.recv()
                try:
                    data = pickle.loads(resp)
                except Exception:
                    # Might be a raw string message (e.g. AUTH_ERROR)
                    try:
                        data = resp.decode('utf-8')
                    except:
                        raise Exception("Invalid response format from server (Pickle truncated/invalid)")

                if isinstance(data, dict) and 'error' in data:
                    raise Exception(data['error'])
                if data == "AUTH_ERROR":
                    raise Exception("Authentication Failed")
                if isinstance(data, str) and data.startswith("ERROR"):
                     raise Exception(data)
                
                if not isinstance(data, dict):
                    raise Exception(f"Server Error: Expected dictionary, received: {data}")

                return data
            else:
                # Reconnect on timeout
                self.socket.close()
                self.socket = self.context.socket(zmq.REQ)
                self.socket.connect(self.server_address if hasattr(self, 'server_address') else self.socket.get(zmq.LAST_ENDPOINT))
                return {'error': 'Communication timeout'}
        except Exception as e:
            return {'error': str(e)}

    def run_variant(self, params):
        cmd = {'command': 'RUN_VARIANT', 'params': params}
        if self.secret: cmd['auth_token'] = self.secret
        try:
            self.socket.send(pickle.dumps(cmd))
            
            poller = zmq.Poller()
            poller.register(self.socket, zmq.POLLIN)
            if poller.poll(30000): # 30s timeout for variant calling
                ret = self.socket.recv()
                try:
                    return pickle.loads(ret)
                except:
                    return []
            else:
                print("DEBUG: Remote variant call timed out")
                return []
        except Exception as e:
            print(f"DEBUG: Remote variant call error: {e}")
            return []

# --- DUPLEX WORKER (On-Demand) ---
class DuplexWorker(QThread):
    progress = pyqtSignal(int)
    results = pyqtSignal(str) # Message
    error = pyqtSignal(str)
    
    def __init__(self, filepath, filters):
        super().__init__()
        self.filepath = filepath
        self.filters = filters # {'min_qs': x, 'min_len': y}
        
    def run(self):
        try:
            import ns_core
            import pandas as pd
            import edlib
            import time
            
            streamer = ns_core.get_streamer(self.filepath)
            
            # 1. Collect Reads (with sequences) that pass filters
            valid_meta = []
            
            min_qs = self.filters.get('min_qs', 0)
            min_len = self.filters.get('min_len', 0)
            
            count = 0
            for batch, meta in streamer.stream_batches(extract_sequences=True):
                for m in meta:
                    if m['qs'] >= min_qs and m['len'] >= min_len:
                        valid_meta.append(m)
                count += len(meta)
                if count % 10000 == 0:
                    self.progress.emit(count)
                    
            if not valid_meta:
                self.results.emit("No reads passed filters.")
                return
                
            df = pd.DataFrame(valid_meta)
            
            # 2. Run Duplex Logic
            has_mapping = False
            if 'rid' in df.columns and 'pos' in df.columns:
                if (df['rid'] != -1).any(): has_mapping = True
                
            def find_duplex_pairs_temporal(group):
                pairs = []
                recs = group.to_dict('records')
                i = 0
                while i < len(recs) - 1:
                    r1 = recs[i]; r2 = recs[i+1]
                    if r2['len'] > r1['len'] + 500: i+=1; continue
                    
                    tail_r1 = r1['tail']; head_r2 = r2['head']
                    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
                    rc_head_r2 = head_r2.translate(trans)[::-1]
                    
                    res_align = edlib.align(tail_r1, rc_head_r2, mode="NW", task="path")
                    if res_align['editDistance'] >= 0:
                        align_len = max(len(tail_r1), len(rc_head_r2))
                        if align_len >= 50:
                            ident = 1.0 - (res_align['editDistance'] / align_len)
                            if ident > 0.70: pairs.append(f"{r1['id']} {r2['id']}"); i+=2; continue
                    i+=1
                return pairs

            def find_duplex_pairs_mapped(group):
                pairs = []
                mapped_recs = [r for r in group.to_dict('records') if r['rid'] != -1]
                mapped_recs.sort(key=lambda x: (x['rid'], x['pos']))
                i = 0
                while i < len(mapped_recs) - 1:
                    r1 = mapped_recs[i]; r2 = mapped_recs[i+1]
                    if r1['rid'] != r2['rid']: i+=1; continue
                    if abs(r1['pos'] - r2['pos']) > 500: i+=1; continue
                    if r1['rev'] == r2['rev']: i+=1; continue
                    
                    res_align = edlib.align(r1['head'], r2['head'], mode="NW", task="path")
                    if res_align['editDistance'] >= 0:
                        align_len = max(len(r1['head']), len(r2['head']))
                        if align_len >= 50:
                            ident = 1.0 - (res_align['editDistance'] / align_len)
                            if ident > 0.70: pairs.append(f"{r1['id']} {r2['id']}"); i+=2; continue
                    i+=1
                return pairs

            all_pairs = []
            if has_mapping:
                results = df.groupby(['ch', 'mx']).apply(find_duplex_pairs_mapped)
            else:
                results = df.groupby(['ch', 'mx']).apply(find_duplex_pairs_temporal)
                
            for p_list in results: all_pairs.extend(p_list)
            
            duplex_pairs = len(all_pairs)
            
            if duplex_pairs > 0:
                fname = f"duplex_candidates_{int(time.time())}.txt"
                with open(fname, "w") as f_out:
                    for p in all_pairs: f_out.write(p + "\n")
                self.results.emit(f"Found {duplex_pairs} pairs. Saved to {fname}")
            else:
                self.results.emit("No duplex pairs found.")
                
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.error.emit(str(e))


class ExportWorker(QThread):
    """Worker thread for exporting reads for selected amplicons."""
    progress = pyqtSignal(int)
    finished = pyqtSignal(int)  # total reads exported
    error = pyqtSignal(str)
    
    def __init__(self, source_file, output_file, amplicon_names, primer_dict, 
                 min_qs, min_len, duplex_only):
        super().__init__()
        self.source_file = source_file
        self.output_file = output_file
        self.amplicon_names = set(amplicon_names)  # Convert to set for faster lookup
        self.primer_dict = primer_dict
        self.min_qs = min_qs
        self.min_len = min_len
        self.duplex_only = duplex_only
    
    def run(self):
        try:
            import pysam
            import ns_amplicon
            
            exported_count = 0
            
            # Determine if source is BAM or FASTQ
            is_bam = self.source_file.endswith('.bam')
            
            # Open source file
            if is_bam:
                source = pysam.AlignmentFile(self.source_file, "rb", check_sq=False)
            else:
                source = pysam.FastxFile(self.source_file)
            
            # Open output file
            if self.output_file.endswith('.bam'):
                if is_bam:
                    output = pysam.AlignmentFile(self.output_file, "wb", template=source)
                else:
                    raise ValueError("Cannot export to BAM from FASTQ source. Please choose FASTQ output.")
            else:
                output = open(self.output_file, 'w')
            
            # Process reads
            total_processed = 0
            
            for read in source:
                total_processed += 1
                
                if total_processed % 10000 == 0:
                    self.progress.emit(total_processed)
                
                # Apply QS filter
                try:
                    if is_bam:
                        qs = read.get_tag("qs") if read.has_tag("qs") else 0
                    else:
                        if read.quality:
                            import numpy as np
                            q_scores = np.frombuffer(read.quality.encode(), dtype=np.int8) - 33
                            qs = np.mean(q_scores)
                        else:
                            qs = 0
                    
                    if qs < self.min_qs:
                        continue
                except:
                    pass
                
                # Apply length filter
                read_len = read.query_length if is_bam else len(read.sequence)
                if read_len < self.min_len:
                    continue
                
                # Apply duplex filter
                if self.duplex_only:
                    try:
                        if is_bam:
                            dx = read.get_tag("dx") if read.has_tag("dx") else 0
                        else:
                            dx = 0
                        if dx != 1:
                            continue
                    except:
                        continue
                
                # Check if read belongs to selected amplicons
                if self.primer_dict:
                    # Primer mode: identify which amplicon this read belongs to
                    amplicon_name = ns_amplicon.identify_primer_for_read(
                        read, self.primer_dict, end_length=150
                    )
                    if amplicon_name and amplicon_name[0] in self.amplicon_names:
                        if self.output_file.endswith('.bam'):
                            output.write(read)
                        else:
                            qual = read.qual if is_bam else read.quality
                            seq = read.query_sequence if is_bam else read.sequence
                            name = read.query_name if is_bam else read.name
                            comment = read.comment if hasattr(read, 'comment') else ""
                            output.write(f"@{name} {comment}\n{seq}\n+\n{qual}\n")
                        exported_count += 1
                else:
                    # Discovery mode: match by position
                    if is_bam and not read.is_unmapped:
                        chrom = read.reference_name
                        start = read.reference_start
                        end = read.reference_end
                        
                        for amp_name in self.amplicon_names:
                            if ":" in amp_name and "-" in amp_name:
                                try:
                                    amp_chrom, coords = amp_name.split(":")
                                    amp_start, amp_end = map(int, coords.split("-"))
                                    
                                    if chrom == amp_chrom and start >= amp_start - 100 and end <= amp_end + 100:
                                        if self.output_file.endswith('.bam'):
                                            output.write(read)
                                        else:
                                            qual = read.qual
                                            seq = read.query_sequence
                                            name = read.query_name
                                            output.write(f"@{name}\n{seq}\n+\n{qual}\n")
                                        exported_count += 1
                                        break
                                except:
                                    pass
            
            source.close()
            if self.output_file.endswith('.bam'):
                output.close()
            else:
                output.close()
            
            self.finished.emit(exported_count)
            
        except Exception as e:
            traceback.print_exc()
            self.error.emit(str(e))

class RemoteDirectoryWatcher(QThread):
    new_files_found = pyqtSignal(list)
    
    def __init__(self, server_address, directory, secret=None, interval=5):
        super().__init__()
        self.server_address = server_address
        self.directory = directory
        self.secret = secret
        self.interval = interval
        self.known_files = set()
        self.running = True
        self.client = RemoteClient(server_address, secret)

    def run(self):
        # Initial scan
        try:
            print(f"RemoteDirectoryWatcher: Scanning {self.directory} (recursive)...")
            items = self.client.list_dir(self.directory, recursive=True)
            print(f"RemoteDirectoryWatcher: Received items type: {type(items)}")
            if isinstance(items, dict) and 'items' in items:
                print(f"RemoteDirectoryWatcher: Found {len(items['items'])} items.")
                initial_files = []
                for item in items['items']:
                   if item['type'] == 'file':
                       initial_files.append(item['path'])
                self.known_files = set(initial_files)
                print(f"RemoteDirectoryWatcher: Found {len(self.known_files)} files initially.")
                if self.known_files:
                    print(f"RemoteDirectoryWatcher: Emitting {len(self.known_files)} files.")
                    self.new_files_found.emit(list(self.known_files))
            else:
                print(f"RemoteDirectoryWatcher: unexpected response: {items}")
        except Exception as e:
            print(f"RemoteDirectoryWatcher Error: {e}")
            import traceback
            traceback.print_exc()
            
        while self.running:
            time.sleep(self.interval)
            try:
                items = self.client.list_dir(self.directory, recursive=True)
                if isinstance(items, dict) and 'items' in items:
                    current_files = set()
                    for item in items['items']:
                        if item['type'] == 'file':
                            current_files.add(item['path'])
                    
                    new_files = list(current_files - self.known_files)
                    if new_files:
                        self.known_files.update(new_files)
                        self.new_files_found.emit(new_files)
            except:
                pass

    def stop(self):
        self.running = False
class RsnapVariantWorker(QThread):
    """
    Worker for running rsnap variant calling and parsing results.
    """
    finished = pyqtSignal(bool, str, list) # success, message, variants
    
    def __init__(self, bam_path, region, af_threshold, output_vcf=None):
        super().__init__()
        self.bam_path = bam_path
        self.region = region # "chr:start-end"
        self.af = af_threshold
        self.output_vcf = output_vcf
        
    def run(self):
        try:
            import subprocess
            
            # 1. Determine Output Filename
            file_to_use = self.output_vcf
            if not file_to_use:
                # Default: chr_start_stop.vcf
                if self.region:
                    try:
                        # expected "chr:start-end"
                        clean = self.region.replace(":", "_").replace("-", "_").replace(",", "")
                        file_to_use = f"{clean}.vcf"
                    except:
                        file_to_use = "temp_variants.vcf"
                else:
                    file_to_use = "temp_variants.vcf"
            
            # 1. Build Command
            # rsnap -b sample.bam -p chr1:1000-2000 --variant --af 0.015 --output [FILE]
            cmd = ["rsnap", "-b", self.bam_path, "--variant", "--af", str(self.af), "--output", file_to_use]
            
            if self.region:
                # Ensure region is clean (no commas)
                clean_region = self.region.replace(",", "")
                cmd.extend(["-p", clean_region])
                
            # 2. Run rsnap
            # Capture output to avoid polluting stdout
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"RsnapVariantWorker rsnap failed: {result.stderr}")
                self.finished.emit(False, f"rsnap failed: {result.stderr}", [])
                return
            else:
                print(f"RsnapVariantWorker finished successfully. Parsing {file_to_use}...")
                
            # 3. Parse VCF
            variants = []
            if os.path.exists(file_to_use):
                with open(file_to_use, "r") as f:
                    for line in f:
                        if line.startswith("#") or not line.strip(): continue
                        
                        # VCF Format: CHROM POS ID REF ALT QUAL FILTER INFO ...
                        parts = line.strip().split("\t")
                        if len(parts) >= 8:
                            chrom = parts[0]
                            try:
                                pos = int(parts[1])
                            except ValueError: continue
                            
                            ref = parts[3]
                            alt = parts[4]
                            info_str = parts[7]
                            
                            # Parse INFO for AF and DP (e.g., "DP=100;AF=0.45;...")
                            af = 0.0
                            dp = 0
                            
                            af = 0.0
                            dp = 0
                            
                            info_dict = {}
                            for item in info_str.split(";"):
                                if "=" in item:
                                    k, v = item.split("=", 1)
                                    info_dict[k] = v
                                    
                            if "AF" in info_dict:
                                try: af = float(info_dict["AF"])
                                except: pass
                            if "DP" in info_dict:
                                try: dp = int(info_dict["DP"])
                                except: pass
                                
                            variants.append({
                                "chrom": chrom,
                                "pos": pos,
                                "ref": ref,
                                "alt": alt,
                                "af": af,
                                "depth": dp,
                                "vaf": af * 100 # For compatibility if needed
                            })
                            
                # Cleanup (Only if user didn't specify output)
                if not self.output_vcf:
                    try:
                        os.remove(file_to_use)
                    except: pass
                
                self.finished.emit(True, f"rsnap found {len(variants)} variants.", variants)
            else:
                 self.finished.emit(False, "No VCF output generated.", [])

        except Exception as e:
            self.finished.emit(False, f"Error running rsnap: {e}", [])

class NanostreamWorker(QThread):
    """
    Worker that calls nanoparse Rust CLI for fast primer matching.
    """
    progress = pyqtSignal(int)
    results = pyqtSignal(object)
    partial_results = pyqtSignal(object)
    finished_file = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, bam_file, primers_path, threads=8, end_length=150, max_edit_dist=3, primer_tolerance=0):
        super().__init__()
        self.bam_file = bam_file
        self.primers_path = primers_path
        self.threads = threads
        self.end_length = end_length
        self.max_edit_dist = max_edit_dist
        self.primer_tolerance = primer_tolerance
        
        # Find unified binary
        self.nanostream_bin = self._find_nanostream()
        # Fallback for legacy
        self.nanostream_bin = self._find_nanoparse() if not self.nanostream_bin else None
        
    def _find_nanostream(self):
        """Locate nanostream binary."""
        import shutil
        found = shutil.which("nanostream")
        if found: return found
        base = os.path.dirname(os.path.abspath(__file__))
        candidates = [
            os.path.join(base, "..", "crates", "nanostream", "target", "release", "nanostream"),
            os.path.join(base, "..", "target", "release", "nanostream"),
        ]
        for c in candidates:
            if os.path.exists(c): return c
        return None

    def _find_nanoparse(self):
        """Locate nanoparse binary."""
        import shutil
        found = shutil.which("nanoparse")
        if found: return found
        base = os.path.dirname(os.path.abspath(__file__))
        candidates = [
            os.path.join(base, "..", "crates", "nanoparse", "target", "release", "nanoparse"),
            os.path.join(base, "..", "target", "release", "nanoparse"),
        ]
        for c in candidates:
            if os.path.exists(c): return c
        return None
        
    def run(self):
        import subprocess
        import json
        
        if not self.nanostream_bin and not self.nanostream_bin:
            self.error.emit("nanostream or nanoparse binary not found. Please build it first.")
            self.finished_file.emit(self.bam_file)
            return
            
        try:
            if self.nanostream_bin:
                cmd = [
                    self.nanostream_bin,
                    "amplicons",
                    self.bam_file,
                    "-p", self.primers_path,
                    "-t", str(self.threads),
                    "--end-length", str(self.end_length),
                    "--max-edit-dist", str(self.max_edit_dist),
                    "--primer-tolerance", str(self.primer_tolerance),
                    "-o", "-"
                ]
            else:
                cmd = [
                    self.nanostream_bin,
                    "amplicons",
                    "-b", self.bam_file,
                    "-p", self.primers_path,
                    "-t", str(self.threads),
                    "--end-length", str(self.end_length),
                    "--max-edit-dist", str(self.max_edit_dist),
                    "--primer-tolerance", str(self.primer_tolerance),
                    "-o", "-"
                ]
            
            self.progress.emit(0)
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode != 0:
                self.error.emit(f"nanoparse failed: {result.stderr}")
                self.finished_file.emit(self.bam_file)
                return
                
            # Parse JSON output
            data = json.loads(result.stdout)
            
            # Convert to nanoMonitor format
            amplicon_stats = {}
            for amp_name, stats in data.get("amplicons", {}).items():
                chrom = stats.get("chrom")
                start = stats.get("start")
                end = stats.get("end")
                region = f"{chrom}:{start}-{end}" if chrom and start is not None and end is not None else None
                
                amplicon_stats[amp_name] = {
                    "count": stats["count"],
                    "median_length": stats["median_length"],
                    "stdev_length": stats.get("std_length", 0),
                    "average_qs": stats.get("avg_qs", 0),
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "region": region,
                    "raw_lengths": [],  # Not needed for display
                }
                
            result_payload = {
                "amplicons": amplicon_stats,  # Changed from amplicon_stats to amplicons
                "unmatched_count": data.get("unmatched_count", 0),
                "total_reads": data.get("total_reads", 0),
                "source": "nanoparse",
            }
            
            self.progress.emit(100)
            self.results.emit(result_payload)
            self.finished_file.emit(self.bam_file)
            
        except subprocess.TimeoutExpired:
            self.error.emit("nanoparse timed out after 10 minutes")
            self.finished_file.emit(self.bam_file)
        except json.JSONDecodeError as e:
            self.error.emit(f"Failed to parse nanoparse output: {e}")
            self.finished_file.emit(self.bam_file)
        except Exception as e:
            self.error.emit(f"NanostreamWorker error: {e}")
            self.finished_file.emit(self.bam_file)
class BatchRsnapWorker(QThread):
    """
    Worker for running rsnap variant calling on multiple amplicons sequentially.
    """
    partial_result = pyqtSignal(str, list) # amplicon_name, variants
    finished = pyqtSignal(bool, str) # success, message
    progress = pyqtSignal(int)
    
    def __init__(self, bam_path, amplicons, af_threshold, genome_path=None, server_address=None, secret=None):
        """
        amplicons: dict {name: region_string}
        """
        super().__init__()
        self.bam_path = bam_path
        self.amplicons = amplicons
        self.af = af_threshold
        self.genome_path = genome_path
        self.server_address = server_address
        self.secret = secret
        self.running = True
        
    def run(self):
        import subprocess
        
        total = len(self.amplicons)
        processed = 0
        
        client = None
        if self.server_address:
            from ns_workers import RemoteClient
            client = RemoteClient(self.server_address, self.secret)

        for name, region in self.amplicons.items():
            if not self.running: break
            
            try:
                # Prepare Temp File
                clean_region = region.replace(":", "_").replace("-", "_").replace(",", "")
                # Create a truly unique temp file to avoid collisions if running repeatedly
                import uuid
                u_id = str(uuid.uuid4())[:8]
                temp_vcf = f"temp_var_{clean_region}_{u_id}.vcf"
                
                if client:
                    # Remote Mode
                    print(f"DEBUG BatchRsnap: Calling remote variant scan on {self.server_address}")
                    variants = client.run_variant({
                        'bam_path': self.bam_path,
                        'region': region,
                        'af': self.af,
                        'genome_path': self.genome_path
                    })
                    print(f"DEBUG BatchRsnap: Received {len(variants)} variants from server")
                else:
                    # Local Mode
                    # Build Command
                    cmd = ["rsnap", "-b", self.bam_path, "-p", region, "--variant", "--af", str(self.af), "--output", temp_vcf]
                    if self.genome_path:
                         cmd.extend(["-r", self.genome_path])  # Use -r for reference, not -g (genes)

                    print(f"DEBUG BatchRsnap: Running cmd: {' '.join(cmd)}")

                    # Run rsnap
                    res = subprocess.run(cmd, capture_output=True, text=True)
                    
                    print(f"DEBUG BatchRsnap: returncode={res.returncode}, stderr={res.stderr[:200] if res.stderr else 'none'}")
                    
                    if res.returncode != 0:
                        print(f"BatchRsnap rsnap failed for {name}: {res.stderr}")
                    
                    variants = []
                    vcf_exists = os.path.exists(temp_vcf)
                    print(f"DEBUG BatchRsnap: VCF exists={vcf_exists}, temp_vcf={temp_vcf}")
                    
                    if res.returncode == 0 and vcf_exists:
                        # Parse VCF
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
                                    
                                    af = 0.0
                                    dp = 0
                                    info_dict = {}
                                    for item in info_str.split(";"):
                                        if "=" in item:
                                            k, v = item.split("=", 1)
                                            info_dict[k] = v
                                    if "AF" in info_dict:
                                        try: af = float(info_dict["AF"])
                                        except: pass
                                    if "DP" in info_dict:
                                        try: dp = int(info_dict["DP"])
                                        except: pass
                                    
                                    variants.append({
                                        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "af": af, "depth": dp, "vaf": af*100
                                    })
                        
                        print(f"DEBUG BatchRsnap: Parsed {len(variants)} variants from VCF")
                        
                        # Clean up
                        try: os.remove(temp_vcf)
                        except: pass
                    
                # Emit result for this amplicon
                self.partial_result.emit(name, variants)
                
            except Exception as e:
                print(f"BatchRsnap error for {name}: {e}")
            
            processed += 1
            self.progress.emit(int(processed / total * 100))
            
        self.finished.emit(True, "Batch Variant Calling Completed")

    def stop(self):
        self.running = False
