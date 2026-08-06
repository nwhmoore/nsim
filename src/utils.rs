//! Shared structure-of-arrays and numerical-accuracy utilities.

/// Three parallel scalar series representing Cartesian components of a vector.
///
/// The vectors are indexed in lockstep. Depending on the owner, an index may
/// identify a particle or a recorded diagnostic sample.
#[derive(Default)]
pub struct VectorSeries {
    /// X components.
    pub x: Vec<f64>,
    /// Y components.
    pub y: Vec<f64>,
    /// Z components.
    pub z: Vec<f64>,
}

/// Accumulates a sum while tracking the rounding error introduced by `f64`
/// addition.
///
/// Compensated summation is useful for diagnostics whose terms can differ
/// greatly in magnitude or partially cancel, such as momentum and potential
/// energy.
#[derive(Debug, Default)]
pub struct KahanAccumulator {
    sum: f64,
    correction: f64,
}

impl KahanAccumulator {
    /// Adds one value to the compensated sum.
    pub fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next_sum = self.sum + adjusted;
        self.correction = (next_sum - self.sum) - adjusted;
        self.sum = next_sum;
    }

    /// Returns the accumulated sum.
    #[must_use]
    pub fn total(&self) -> f64 {
        self.sum
    }
}
