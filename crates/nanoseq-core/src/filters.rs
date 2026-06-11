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

pub fn parse_duration(s: &str) -> Result<chrono::Duration> {
    let s = s.trim();
    if s.is_empty() {
        return Err(anyhow!("Empty duration string"));
    }

    let is_negative = s.starts_with('-');
    let abs_s = if is_negative { &s[1..] } else { s };

    let mut total_seconds = 0i64;
    let mut current_val = 0i64;
    let mut has_digit = false;
    let mut has_unit = false;

    if abs_s.chars().all(|c| c.is_ascii_digit()) {
        let val = abs_s
            .parse::<i64>()
            .map_err(|_| anyhow!("Invalid number: {abs_s}"))?;
        total_seconds = val * 60; // default to minutes
    } else {
        let mut chars = abs_s.chars().peekable();
        while let Some(c) = chars.next() {
            if c.is_ascii_digit() {
                current_val = current_val * 10 + c.to_digit(10).unwrap() as i64;
                has_digit = true;
            } else if c.is_alphabetic() {
                if !has_digit {
                    return Err(anyhow!("Missing number before unit '{c}'"));
                }
                let unit = c.to_ascii_lowercase();
                match unit {
                    'd' => total_seconds += current_val * 24 * 3600,
                    'h' => total_seconds += current_val * 3600,
                    'm' => total_seconds += current_val * 60,
                    's' => total_seconds += current_val,
                    _ => return Err(anyhow!("Unknown unit '{unit}'")),
                }
                current_val = 0;
                has_digit = false;
                has_unit = true;
            } else if c.is_whitespace() {
                continue;
            } else {
                return Err(anyhow!("Invalid character in duration: '{c}'"));
            }
        }
        if has_digit && !has_unit {
            total_seconds += current_val * 60; // fallback to minutes
        }
    }

    let seconds = if is_negative { -total_seconds } else { total_seconds };
    Ok(chrono::Duration::seconds(seconds))
}

pub fn resolve_relative_window(
    time_start: &str,
    time_end: &str,
    run_start: DateTime<FixedOffset>,
    run_end: DateTime<FixedOffset>,
) -> Result<TimeWindow> {
    // 1. Try to parse time_start as absolute, or parse it as a relative duration
    let start = if let Some(dt) = DateTime::parse_from_rfc3339(time_start.trim()).ok() {
        dt
    } else {
        let dur = parse_duration(time_start)?;
        if dur.num_seconds() >= 0 {
            run_start + dur
        } else {
            run_end + dur
        }
    };

    // 2. Try to parse time_end as absolute, or parse it as a relative duration
    let end = if let Some(dt) = DateTime::parse_from_rfc3339(time_end.trim()).ok() {
        dt
    } else {
        let dur = parse_duration(time_end)?;
        if dur.num_seconds() >= 0 {
            run_start + dur
        } else {
            run_end + dur
        }
    };

    if start > end {
        return Err(anyhow!(
            "Resolved time window is invalid: start ({}) is after end ({})",
            start.to_rfc3339(),
            end.to_rfc3339()
        ));
    }

    Ok(TimeWindow { start, end })
}

#[cfg(test)]
mod tests {
    use super::{parse_pore_range, parse_time_window, parse_duration, resolve_relative_window};
    use chrono::DateTime;

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

    #[test]
    fn parses_durations() {
        assert_eq!(parse_duration("0").unwrap().num_seconds(), 0);
        assert_eq!(parse_duration("90").unwrap().num_seconds(), 5400);
        assert_eq!(parse_duration("2h30m").unwrap().num_seconds(), 9000);
        assert_eq!(parse_duration("10s").unwrap().num_seconds(), 10);
        assert_eq!(parse_duration("-1h").unwrap().num_seconds(), -3600);
    }

    #[test]
    fn resolves_relative_windows() {
        let start = DateTime::parse_from_rfc3339("2026-03-17T12:00:00+00:00").unwrap();
        let end = DateTime::parse_from_rfc3339("2026-03-17T15:00:00+00:00").unwrap();

        // 1. Durations relative to start/end
        let window1 = resolve_relative_window("0", "90m", start, end).unwrap();
        assert_eq!(window1.start, start);
        assert_eq!(window1.end, start + chrono::Duration::minutes(90));

        // 2. Negative offset relative to end
        let window2 = resolve_relative_window("1h", "-1h", start, end).unwrap();
        assert_eq!(window2.start, start + chrono::Duration::hours(1));
        assert_eq!(window2.end, end - chrono::Duration::hours(1));

        // 3. Absolute start, relative end
        let window3 = resolve_relative_window("2026-03-17T12:30:00+00:00", "120m", start, end).unwrap();
        assert_eq!(window3.start, DateTime::parse_from_rfc3339("2026-03-17T12:30:00+00:00").unwrap());
        assert_eq!(window3.end, start + chrono::Duration::minutes(120));
    }
}
