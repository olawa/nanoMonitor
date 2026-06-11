#!/bin/bash
# merge_and_filter.sh
# Merges raw barcode FASTQ files from a nanopore run and runs nanostream amplicons filtering.

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -s, --samples PATH       Path to samples.tsv (default: ./samples.tsv)"
    echo "  -f, --fastq-pass PATH    Path to fastq_pass directory (default: ./fastq_pass)"
    echo "  -p, --primers PATH       Path to primers.tsv (default: ./primers.tsv)"
    echo "  -l, --len RANGE          Length filter (e.g. 1000-3000 or 1000) (default: 0)"
    echo "  -q, --min-qs SCORE       Minimum mean quality score filter (default: 0)"
    echo "  -o, --output-dir PATH    Output directory (default: ./output_filtered)"
    echo "  -t, --threads NUM        Number of processing threads (default: 8)"
    echo "  --split                  Split matched reads by amplicon name"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "The samples.tsv format should be tab-separated:"
    echo "  run_id <tab> sample_id <tab> barcode"
}

# Default values
SAMPLES="./samples.tsv"
FASTQ_PASS="./fastq_pass"
PRIMERS="./primers.tsv"
LEN="0"
MIN_QS="0"
OUTPUT_DIR="./output_filtered"
THREADS="8"
SPLIT_BY_AMPLICON=false

# Parse command line options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samples)
            SAMPLES="$2"
            shift 2
            ;;
        -f|--fastq-pass)
            FASTQ_PASS="$2"
            shift 2
            ;;
        -p|--primers)
            PRIMERS="$2"
            shift 2
            ;;
        -l|--len)
            LEN="$2"
            shift 2
            ;;
        -q|--min-qs)
            MIN_QS="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --split)
            SPLIT_BY_AMPLICON=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# Ensure nanostream binary is built
if [ ! -f "./target/release/nanostream" ]; then
    echo -e "${YELLOW}nanostream release binary not found. Building now...${NC}"
    cargo build --release --bin nanostream
fi

# Check paths
if [ ! -f "$SAMPLES" ]; then
    echo -e "${RED}Error: Samples file not found at $SAMPLES${NC}"
    exit 1
fi

if [ ! -d "$FASTQ_PASS" ]; then
    echo -e "${RED}Error: fastq_pass directory not found at $FASTQ_PASS${NC}"
    exit 1
fi

if [ ! -f "$PRIMERS" ]; then
    echo -e "${RED}Error: Primers file not found at $PRIMERS${NC}"
    exit 1
fi

# Create directories
MERGED_DIR="$OUTPUT_DIR/merged"
mkdir -p "$MERGED_DIR"

echo -e "${BLUE}=== Starting Merge & Filter Pipeline ===${NC}"
echo "Samples sheet:  $SAMPLES"
echo "fastq_pass:     $FASTQ_PASS"
echo "Primers:        $PRIMERS"
echo "Length filter:  $LEN"
echo "Quality filter: $MIN_QS"
echo "Output dir:     $OUTPUT_DIR"
echo "Threads:        $THREADS"
echo "Split output:   $SPLIT_BY_AMPLICON"
echo "========================================"

# Process samples.tsv line by line
# Skip empty lines and comments (lines starting with #)
line_num=0
while IFS=$'\t' read -r run_id sample_id barcode || [ -n "$run_id" ]; do
    # Trim whitespace
    run_id=$(echo "$run_id" | xargs)
    sample_id=$(echo "$sample_id" | xargs)
    barcode=$(echo "$barcode" | xargs)

    # Skip header, empty lines or comments
    if [[ -z "$run_id" || "$run_id" =~ ^# || "$run_id" == "run_id" ]]; then
        continue
    fi

    line_num=$((line_num + 1))
    echo -e "\n${YELLOW}[$line_num] Processing sample: $sample_id (Run: $run_id, Barcode: $barcode)${NC}"

    # Find matching FASTQ files
    FILES=()
    
    # 1. Check run_id/barcode subfolder
    if [ -d "$FASTQ_PASS/$run_id/$barcode" ]; then
        while IFS= read -r -d '' file; do
            FILES+=("$file")
        done < <(find "$FASTQ_PASS/$run_id/$barcode" -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) -print0)
    # 2. Check barcode subfolder directly
    elif [ -d "$FASTQ_PASS/$barcode" ]; then
        while IFS= read -r -d '' file; do
            FILES+=("$file")
        done < <(find "$FASTQ_PASS/$barcode" -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) -print0)
    # 3. Scan directory for matching names
    else
        while IFS= read -r -d '' file; do
            FILES+=("$file")
        done < <(find "$FASTQ_PASS" -type f \( -name "*${barcode}*.fastq" -o -name "*${barcode}*.fastq.gz" -o -name "*${barcode}*.fq" -o -name "*${barcode}*.fq.gz" \) -print0)
    fi

    if [ ${#FILES[@]} -eq 0 ]; then
        echo -e "${RED}[ERROR] No FASTQ files found for run $run_id, barcode $barcode (sample $sample_id). Skipping.${NC}"
        continue
    fi

    MERGED_FILE="$MERGED_DIR/${sample_id}.fastq.gz"
    echo -e "${GREEN}[MERGE]${NC} Merging ${#FILES[@]} file(s) into $MERGED_FILE..."
    
    # Check if first file is compressed
    first_file="${FILES[0]}"
    if [[ "$first_file" =~ \.gz$ ]]; then
        cat "${FILES[@]}" > "$MERGED_FILE"
    else
        cat "${FILES[@]}" | gzip -c > "$MERGED_FILE"
    fi

    # Run nanostream amplicons
    CMD=("./target/release/nanostream" "amplicons" \
         "--primers" "$PRIMERS" \
         "--len" "$LEN" \
         "--min-qs" "$MIN_QS" \
         "--threads" "$THREADS" \
         "--output" "$OUTPUT_DIR/${sample_id}_stats.json" \
         "--output-fastq" "$OUTPUT_DIR/${sample_id}_matched.fastq.gz" \
         "--output-dimers" "$OUTPUT_DIR/${sample_id}_dimers.fastq.gz" \
         "--summary")

    if [ "$SPLIT_BY_AMPLICON" = true ]; then
        CMD+=("--split-by-amplicon")
    fi

    CMD+=("$MERGED_FILE")

    echo -e "${GREEN}[RUN]${NC} Running nanostream amplicons..."
    if "${CMD[@]}"; then
        echo -e "${GREEN}[SUCCESS]${NC} Sample $sample_id processed successfully."
    else
        echo -e "${RED}[ERROR]${NC} Failed to process sample $sample_id."
    fi
done < "$SAMPLES"

echo -e "\n${BLUE}=== Pipeline Completed ===${NC}"
