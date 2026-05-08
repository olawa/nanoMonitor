use pyo3.prelude::*;
use std::collections::BinaryHeap;
use std::cmp::Reverse;

#[pyfunction]
fn calculate_lanes(reads: Vec<(u32, u32)>) -> Vec<usize> {
    // Min-heap storing (end_position, lane_index)
    // We use Reverse because Rust's BinaryHeap is a Max-Heap
    let mut heap: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
    let mut lanes = Vec::with_capacity(reads.len());
    let mut next_lane = 0;

    // Reads must be sorted by start position before calling this!
    // (We assume Python does the sort as it's fast enough there)

    for (start, end) in reads {
        let mut assigned_lane = 0;

        // Check if the earliest ending lane is free
        if let Some(Reverse((lane_end, lane_idx))) = heap.peek() {
            if *lane_end < start {
                assigned_lane = *lane_idx;
                heap.pop(); // Remove from heap to update it
            } else {
                assigned_lane = next_lane;
                next_lane += 1;
            }
        } else {
            assigned_lane = next_lane;
            next_lane += 1;
        }

        lanes.push(assigned_lane);
        // Add back with new end position + gap (e.g., 5bp)
        heap.push(Reverse((end + 5, assigned_lane)));
    }

    lanes
}

#[pymodule]
fn ns_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calculate_lanes, m)?)?;
    Ok(())
}