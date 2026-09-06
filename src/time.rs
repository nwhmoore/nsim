//! Internal clock, timestep, and diagnostic scheduling.

/// Simulation clock, timestep, and diagnostic schedule.
#[derive(Clone)]
pub struct Time {
    /// Current simulation time.
    pub(crate) current: f64,
    /// Simulation timestep.
    pub(crate) step: f64,
    /// Diagnostic recording schedule.
    pub(crate) diagnostic_schedule: DiagnosticSchedule,
}

impl Default for Time {
    fn default() -> Self {
        Self {
            current: 0.0,
            step: 0.001,
            diagnostic_schedule: DiagnosticSchedule::default(),
        }
    }
}

impl Time {
    pub(crate) fn set_diagnostic_interval(&mut self, dt: f64) {
        self.diagnostic_schedule.diagnostic_interval = dt;
        self.diagnostic_schedule.next_diagnostic_record =
            self.diagnostic_schedule.diagnostic_interval;
    }
}

#[derive(Clone)]
pub(crate) struct DiagnosticSchedule {
    pub(crate) diagnostic_interval: f64,
    pub(crate) next_diagnostic_record: f64,
}

impl Default for DiagnosticSchedule {
    fn default() -> Self {
        DiagnosticSchedule {
            diagnostic_interval: f64::INFINITY,
            next_diagnostic_record: f64::INFINITY,
        }
    }
}
