const RC_TABLE: [u8; 256] = {
    let mut table = [b'N'; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = match i as u8 {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        };
        i += 1;
    }
    table
};

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| RC_TABLE[b as usize]).collect()
}

pub fn hamming_distance(a: &[u8], b: &[u8]) -> usize {
    if a.len() != b.len() {
        return usize::MAX;
    }
    a.iter().zip(b.iter()).filter(|&(x, y)| x != y).count()
}

pub fn edit_distance(a: &[u8], b: &[u8], max_dist: usize) -> usize {
    let n = a.len();
    let m = b.len();

    if n.abs_diff(m) > max_dist {
        return usize::MAX;
    }
    if n == 0 {
        return m;
    }
    if m == 0 {
        return n;
    }
    
    // Simple iterative Levenshtein with early exit
    let mut prev = (0..=m).collect::<Vec<_>>();
    let mut curr = vec![0usize; m + 1];

    for i in 1..=n {
        curr[0] = i;
        let mut min_in_row = i;
        for j in 1..=m {
            let cost = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            let val = (prev[j] + 1).min(curr[j - 1] + 1).min(prev[j - 1] + cost);
            curr[j] = val;
            min_in_row = min_in_row.min(val);
        }
        if min_in_row > max_dist {
            return usize::MAX;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let final_dist = prev[m];
    if final_dist <= max_dist {
        final_dist
    } else {
        usize::MAX
    }
}
