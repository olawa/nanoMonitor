#!/usr/bin/env python3
"""
Test script for parse_amplicon_name function.
Tests various amplicon name formats to ensure correct parsing.
"""

import sys
sys.path.insert(0, '/Users/olwal516/dev/Genomics_Suite/apps_python/nanoStream')

from ns_amplicon import parse_amplicon_name

def test_parse_amplicon_name():
    """Test the parse_amplicon_name function with various inputs."""
    
    test_cases = [
        {
            "input": "chr17:41196312-41277500(BRCA1_ex2-5)",
            "expected": {
                "chrom": "chr17",
                "start": 41196312,
                "end": 41277500,
                "gene_name": "BRCA1_ex2-5",
                "genes": ["BRCA1"]
            }
        },
        {
            "input": "chr17:41196312-41277500(BRCA1_ex2-5,TP53_ex1)",
            "expected": {
                "chrom": "chr17",
                "start": 41196312,
                "end": 41277500,
                "gene_name": "BRCA1_ex2-5,TP53_ex1",
                "genes": ["BRCA1", "TP53"]
            }
        },
        {
            "input": "chr17:41196312-41277500",
            "expected": {
                "chrom": "chr17",
                "start": 41196312,
                "end": 41277500,
                "gene_name": None,
                "genes": []
            }
        },
        {
            "input": "PRIMER1...PRIMER2",
            "expected": {
                "chrom": None,
                "start": None,
                "end": None,
                "gene_name": None,
                "genes": []
            }
        },
        {
            "input": "13:32889611-32973805(BRCA2_ex10-11)",
            "expected": {
                "chrom": "13",
                "start": 32889611,
                "end": 32973805,
                "gene_name": "BRCA2_ex10-11",
                "genes": ["BRCA2"]
            }
        }
    ]
    
    print("Testing parse_amplicon_name function...\n")
    
    passed = 0
    failed = 0
    
    for i, test in enumerate(test_cases, 1):
        input_name = test["input"]
        expected = test["expected"]
        result = parse_amplicon_name(input_name)
        
        # Check if result matches expected
        match = True
        errors = []
        
        for key in expected:
            if result[key] != expected[key]:
                match = False
                errors.append(f"  {key}: expected {expected[key]}, got {result[key]}")
        
        if match:
            print(f"✓ Test {i}: PASSED")
            print(f"  Input: {input_name}")
            passed += 1
        else:
            print(f"✗ Test {i}: FAILED")
            print(f"  Input: {input_name}")
            for error in errors:
                print(error)
            failed += 1
        print()
    
    print(f"\nResults: {passed} passed, {failed} failed out of {len(test_cases)} tests")
    
    return failed == 0

if __name__ == "__main__":
    success = test_parse_amplicon_name()
    sys.exit(0 if success else 1)
