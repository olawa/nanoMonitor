use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RemoteStatus {
    Disconnected,
    Connecting,
    Connected,
}

impl RemoteStatus {
    pub fn label(self) -> &'static str {
        match self {
            RemoteStatus::Disconnected => "Disconnected",
            RemoteStatus::Connecting => "Connecting",
            RemoteStatus::Connected => "Connected",
        }
    }
}

#[derive(Debug, Clone)]
pub struct RemoteConfig {
    pub enabled: bool,
    pub endpoint: String,
    pub auth_token: String,
    pub status: RemoteStatus,
}

impl Default for RemoteConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            endpoint: "tcp://127.0.0.1:5555".into(),
            auth_token: String::new(),
            status: RemoteStatus::Disconnected,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum MonitorRequest {
    Ping,
    Start {
        input_path: String,
        mode: String,
        max_reads: u32,
    },
    Stop,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum MonitorEvent {
    Pong,
    Progress {
        reads_processed: u64,
        percent: f32,
    },
    ResultSummary {
        total_reads: u64,
        filtered_reads: u64,
    },
    Error {
        message: String,
    },
}
