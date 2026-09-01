//! contains the Kahan accumulator to account for floating point error

use crate::math_util::Vector3;
use std::f64;

/// Accumulates a sum while tracking the rounding error introduced by `f64`
/// addition.
///
/// Compensated summation is useful when terms can differ greatly in magnitude
/// or partially cancel, such as force components, momentum, and potential
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

    /// Resets the accumulated total and compensation error while retaining
    /// the accumulator's storage.
    pub fn reset(&mut self) {
        self.sum = 0.0;
        self.correction = 0.0;
    }
}

/// Three parallel accumulators representing Cartesian components of a vector.
#[derive(Default)]
pub struct Kahan3 {
    x: KahanAccumulator,
    y: KahanAccumulator,
    z: KahanAccumulator,
}

impl Kahan3 {
    /// Adds one vector-valued sample to the compensated total in each axis.
    pub fn add(&mut self, x: f64, y: f64, z: f64) {
        self.x.add(x);
        self.y.add(y);
        self.z.add(z);
    }

    /// Returns the compensated total for the X, Y, and Z components.
    #[must_use]
    pub fn total(&self) -> Vector3 {
        Vector3 {
            x: self.x.total(),
            y: self.y.total(),
            z: self.z.total(),
        }
    }
}
