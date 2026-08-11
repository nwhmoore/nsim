//! Shared structure-of-arrays and numerical-accuracy utilities.

/// Three Cartesian components of a vector.
pub struct Vector3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

/// Pairwise geometric data.
pub struct Geometry {
    /// Relative displacement vector between two particles.
    pub r_vec: Vector3,
    /// Euclidean separation magnitude.
    pub dist: f64,
    /// Inverse cube of the separation magnitude.
    pub inv_dist_cubed: f64,
}

impl Geometry {
    /// Computes the relative geometry for one particle pair.
    pub fn calculate_geometry(r_vec: Vector3) -> Self {
        let dist_squared = r_vec.x * r_vec.x + r_vec.y * r_vec.y + r_vec.z * r_vec.z;
        let dist = dist_squared.sqrt();
        let inv_dist_cubed = 1.0 / (dist_squared * dist);

        Geometry {
            r_vec,
            dist,
            inv_dist_cubed,
        }
    }
}

/// Three parallel scalar series representing Cartesian components of a vector.
///
/// The vectors are indexed in lockstep. Depending on the owner, an index may
/// identify a particle or a recorded diagnostic sample.
#[derive(Default)]
pub struct Vector3Series {
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
