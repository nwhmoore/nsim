//! `nsim` is a small Rust N-body simulation prototype. It features a modular
//! integrator and force system.

#![warn(missing_docs)]
#![warn(clippy::all, clippy::pedantic)]
#![warn(rustdoc::broken_intra_doc_links)]

pub mod diagnostics;
pub mod force;
pub mod integration;
pub mod math_util;
pub mod particle;
pub mod simulation;
pub mod time;
