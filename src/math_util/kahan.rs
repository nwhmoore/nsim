use crate::math_util::Vector3;

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

#[derive(Default)]
pub struct Kahan3 {
    x: KahanAccumulator,
    y: KahanAccumulator,
    z: KahanAccumulator,
}

impl Kahan3 {
    /// Adds one vector-valued sample to the compensated total in each axis.
    pub fn add(&mut self, value: &Vector3) {
        self.x.add(value.x);
        self.y.add(value.y);
        self.z.add(value.z);
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

/// Per-particle Kahan-compensated acceleration totals used during one force
/// evaluation.
pub struct Kahan3Series {
    /// Accumulated X-component accelerations for each particle.
    pub x: Vec<KahanAccumulator>,
    /// Accumulated Y-component accelerations for each particle.
    pub y: Vec<KahanAccumulator>,
    /// Accumulated Z-component accelerations for each particle.
    pub z: Vec<KahanAccumulator>,
}

impl Kahan3Series {
    /// Creates an accumulator with one compensated total per `capacity`.
    pub fn with_len(capacity: usize) -> Self {
        Self {
            x: (0..capacity).map(|_| KahanAccumulator::default()).collect(),
            y: (0..capacity).map(|_| KahanAccumulator::default()).collect(),
            z: (0..capacity).map(|_| KahanAccumulator::default()).collect(),
        }
    }

    pub fn reset_at(&mut self, idx: usize) {
        self.x[idx].reset();
        self.y[idx].reset();
        self.z[idx].reset();
    }

    pub fn len(&self) -> usize {
        debug_assert_eq!(self.x.len(), self.y.len());
        debug_assert_eq!(self.x.len(), self.z.len());

        self.x.len()
    }

    #[allow(unused)]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}
