use crate::math_util::vector3::Vector3;

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
    pub fn add(&mut self, value: &Vector3) {
        self.x.add(value.x);
        self.y.add(value.y);
        self.z.add(value.z);
    }

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
    x: Vec<KahanAccumulator>,
    /// Accumulated Y-component accelerations for each particle.
    y: Vec<KahanAccumulator>,
    /// Accumulated Z-component accelerations for each particle.
    z: Vec<KahanAccumulator>,
}

impl Kahan3Series {
    /// Creates an accumulator with one compensated total per particle.
    pub fn new(number_particles: usize) -> Self {
        Self {
            x: (0..number_particles)
                .map(|_| KahanAccumulator::default())
                .collect(),
            y: (0..number_particles)
                .map(|_| KahanAccumulator::default())
                .collect(),
            z: (0..number_particles)
                .map(|_| KahanAccumulator::default())
                .collect(),
        }
    }

    /// Adds one acceleration contribution to the stored total for a particle.
    pub fn add(&mut self, particle_idx: usize, acceleration: &Vector3) {
        self.x[particle_idx].add(acceleration.x);
        self.y[particle_idx].add(acceleration.y);
        self.z[particle_idx].add(acceleration.z);
    }

    pub fn total(&self, particle_idx: usize) -> Vector3 {
        Vector3 {
            x: self.x[particle_idx].total(),
            y: self.y[particle_idx].total(),
            z: self.z[particle_idx].total(),
        }
    }

    pub fn reset_at(&mut self, idx: usize) {
        self.x[idx].reset();
        self.y[idx].reset();
        self.z[idx].reset();
    }
}
