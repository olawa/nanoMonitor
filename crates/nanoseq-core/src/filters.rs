use anyhow::{anyhow, Result};
use chrono::{DateTime, FixedOffset};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PoreRange {
    pub min: i64,
    pub max: i64,
}

impl PoreRange {
    pub fn contains(&self, value: i64) -> bool {
        value >= self.min && value <= self.max
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TimeWindow {
    pub start: DateTime<FixedOffset>,
    pub end: DateTime<FixedOffset>,
}

impl TimeWindow {
    pub fn contains(&self, value: &DateTime<FixedOffset>) -> bool {
        value >= &self.start && value <= &self.end
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ReadFilter {
    pub min_qv: Option<f64>,
    pub min_len: Option<usize>,
    pub max_len: Option<usize>,
    pub duplex_only: bool,
    pub pore_range: Option<PoreRange>,
    pub time_window: Option<TimeWindow>,
    pub max_reads: Option<usize>,
}

pub fn parse_pore_range(raw: &str) -> Result<PoreRange> {
    let trimmed = raw.trim();
    let (start, end) = trimmed
        .split_once('-')
        .ok_or_else(|| anyhow!("Invalid range '{trimmed}', expected START-END"))?;
    let min = start
        .trim()
        .parse::<i64>()
        .map_err(|_| anyhow!("Invalid range start in '{trimmed}'"))?;
    let max = end
        .trim()
        .parse::<i64>()
        .map_err(|_| anyhow!("Invalid range end in '{trimmed}'"))?;
    if min > max {
        return Err(anyhow!("Invalid range '{trimmed}': start > end"));
    }
    Ok(PoreRange { min, max })
}

pub fn parse_time_window(start: &str, end: &str) -> Result<TimeWindow> {
    let start = DateTime::parse_from_rfc3339(start.trim())
        .map_err(|e| anyhow!("Invalid start time: {e}"))?;
    let end =
        DateTime::parse_from_rfc3339(end.trim()).map_err(|e| anyhow!("Invalid end time: {e}"))?;
    if start > end {
        return Err(anyhow!("Invalid time window: start > end"));
    }
    Ok(TimeWindow { start, end })
}

#[cfg(test)]
mod tests {
    use super::{parse_pore_range, parse_time_window};

    #[test]
    fn parses_pore_range() {
        let range = parse_pore_range("1-2500").unwrap();
        assert_eq!(range.min, 1);
        assert_eq!(range.max, 2500);
        assert!(range.contains(99));
        assert!(!range.contains(3000));
    }

    #[test]
    fn parses_time_window() {
        let window =
            parse_time_window("2026-03-17T10:00:00+00:00", "2026-03-17T11:00:00+00:00").unwrap();
        assert!(window.start < window.end);
    }
}
