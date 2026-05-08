import time
import ns_amplicon
import ns_core

def test_parallel_vs_sequential():
    """
    Test parallel processing vs sequential.
    """
    # Mock primers
    primers = {
        "TEST_FWD": "ACGTACGTACGTACGTACGT",
        "TEST_REV": "TGCATGCATGCATGCATGCA"
    }
    
    # Test file (use a real BAM/FASTQ if available)
    test_file = "test.bam"  # Replace with actual file
    
    print("Testing Parallel Processing...")
    print("=" * 50)
    
    # Sequential (threads=1)
    print("\n1. Sequential (threads=1):")
    start = time.time()
    result_seq = ns_amplicon.run_analysis(
        test_file, primers, None, {}, threads=1, 
        progress_cb=None, partial_cb=None
    )
    time_seq = time.time() - start
    print(f"   Time: {time_seq:.2f}s")
    print(f"   Reads: {result_seq['summary']['total_reads_processed']}")
    
    # Parallel (threads=8)
    print("\n2. Parallel (threads=8):")
    start = time.time()
    result_par = ns_amplicon.run_analysis(
        test_file, primers, None, {}, threads=8, 
        progress_cb=None, partial_cb=None
    )
    time_par = time.time() - start
    print(f"   Time: {time_par:.2f}s")
    print(f"   Reads: {result_par['summary']['total_reads_processed']}")
    
    # Speedup
    speedup = time_seq / time_par
    print(f"\n3. Speedup: {speedup:.2f}x")
    
    # Verify results are identical
    print("\n4. Verifying results...")
    seq_amps = result_seq.get('amplicons', {})
    par_amps = result_par.get('amplicons', {})
    
    if seq_amps.keys() == par_amps.keys():
        print("   ✓ Amplicon names match")
        for name in seq_amps:
            if seq_amps[name]['count'] == par_amps[name]['count']:
                print(f"   ✓ {name}: counts match ({seq_amps[name]['count']})")
            else:
                print(f"   ✗ {name}: counts differ! Seq={seq_amps[name]['count']}, Par={par_amps[name]['count']}")
    else:
        print("   ✗ Amplicon names differ!")
        print(f"   Sequential: {seq_amps.keys()}")
        print(f"   Parallel: {par_amps.keys()}")
    
    print("\n" + "=" * 50)
    print("Test Complete!")

if __name__ == "__main__":
    # Note: This test requires a real BAM/FASTQ file
    # Replace "test.bam" with an actual file path
    print("Note: Update test_file path before running.")
    # test_parallel_vs_sequential()
