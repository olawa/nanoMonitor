use serde::{Deserialize, Serialize};
use std::io::{BufRead, BufReader, Write};
use std::net::{TcpStream, ToSocketAddrs};
use std::sync::mpsc::{channel, Sender};
use std::thread;
use std::time::Duration;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RemoteStatus {
    Disconnected,
    Connecting,
    Connected,
}

impl RemoteStatus {
    pub fn label(self) -> &'static str {
        match self {
            RemoteStatus::Disconnected => "Offline",
            RemoteStatus::Connecting => "Connecting...",
            RemoteStatus::Connected => "Connected (LAN)",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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
            endpoint: "127.0.0.1:5555".into(),
            auth_token: String::new(),
            status: RemoteStatus::Disconnected,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
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

#[derive(Debug, Serialize, Deserialize, Clone)]
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

/// Spawns a background thread to handle TCP connection to the remote LAN node.
/// Returns a control channel sender to push requests to the remote node, and a receiver to disconnect.
pub fn spawn_remote_client(
    endpoint: String,
    auth_token: String,
    event_tx: Sender<MonitorEvent>,
) -> (Sender<MonitorRequest>, Sender<()>) {
    let (req_tx, req_rx) = channel::<MonitorRequest>();
    let (stop_tx, stop_rx) = channel::<()>();

    thread::spawn(move || {
        // Parse endpoint address (strip tcp:// if present)
        let addr_str = endpoint.trim().trim_start_matches("tcp://");
        let socket_addr = match addr_str.to_socket_addrs() {
            Ok(mut addrs) => match addrs.next() {
                Some(addr) => addr,
                None => {
                    let _ = event_tx.send(MonitorEvent::Error {
                        message: "Invalid remote address: no resolved addresses".into(),
                    });
                    return;
                }
            },
            Err(e) => {
                let _ = event_tx.send(MonitorEvent::Error {
                    message: format!("Invalid remote address: {}", e),
                });
                return;
            }
        };

        // Attempt connection with timeout
        let stream = match TcpStream::connect_timeout(&socket_addr, Duration::from_secs(4)) {
            Ok(s) => s,
            Err(e) => {
                let _ = event_tx.send(MonitorEvent::Error {
                    message: format!("Connection failed: {}", e),
                });
                return;
            }
        };

        // Sockets successfully connected! Configure socket timeouts
        let _ = stream.set_read_timeout(Some(Duration::from_millis(100)));
        let _ = stream.set_write_timeout(Some(Duration::from_secs(2)));

        let read_stream = match stream.try_clone() {
            Ok(s) => s,
            Err(e) => {
                let _ = event_tx.send(MonitorEvent::Error {
                    message: format!("Socket cloning failed: {}", e),
                });
                return;
            }
        };

        let mut write_stream = stream;

        // Perform authentication handshake if a token is supplied
        if !auth_token.trim().is_empty() {
            let handshake_req = format!("AUTH {}\n", auth_token.trim());
            if let Err(e) = write_stream.write_all(handshake_req.as_bytes()) {
                let _ = event_tx.send(MonitorEvent::Error {
                    message: format!("Handshake write failed: {}", e),
                });
                return;
            }
        }

        // Successfully connected! Notify UI
        let _ = event_tx.send(MonitorEvent::Pong);

        // Read thread loop
        let reader_event_tx = event_tx.clone();
        let read_handle = thread::spawn(move || {
            let mut reader = BufReader::new(read_stream);
            let mut line = String::new();
            loop {
                line.clear();
                match reader.read_line(&mut line) {
                    Ok(0) => {
                        // Connection closed by server
                        let _ = reader_event_tx.send(MonitorEvent::Error {
                            message: "Connection closed by remote node".into(),
                        });
                        break;
                    }
                    Ok(_) => {
                        let trimmed = line.trim();
                        if trimmed.is_empty() {
                            continue;
                        }
                        if let Ok(event) = serde_json::from_str::<MonitorEvent>(trimmed) {
                            if reader_event_tx.send(event).is_err() {
                                break; // Receiver hung up
                            }
                        } else {
                            // Non-JSON logging fallback or raw message
                            let _ = reader_event_tx.send(MonitorEvent::Error {
                                message: format!("Raw message: {}", trimmed),
                            });
                        }
                    }
                    Err(ref e) if e.kind() == std::io::ErrorKind::WouldBlock || e.kind() == std::io::ErrorKind::TimedOut => {
                        // Check again, timeout occurred
                        thread::sleep(Duration::from_millis(50));
                    }
                    Err(e) => {
                        let _ = reader_event_tx.send(MonitorEvent::Error {
                            message: format!("Read error: {}", e),
                        });
                        break;
                    }
                }
            }
        });

        // Write thread / control loop
        loop {
            // Check stop signal
            if stop_rx.try_recv().is_ok() {
                break;
            }

            // Check if there are requests to send
            if let Ok(req) = req_rx.try_recv() {
                if let Ok(mut serialized) = serde_json::to_string(&req) {
                    serialized.push('\n');
                    if let Err(e) = write_stream.write_all(serialized.as_bytes()) {
                        let _ = event_tx.send(MonitorEvent::Error {
                            message: format!("Write error: {}", e),
                        });
                        break;
                    }
                    let _ = write_stream.flush();
                }
            }

            thread::sleep(Duration::from_millis(50));
        }

        // Clean up connections
        let _ = write_stream.shutdown(std::net::Shutdown::Both);
        let _ = read_handle.join();
    });

    (req_tx, stop_tx)
}
