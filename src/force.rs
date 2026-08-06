//! Gravitational acceleration calculation and reusable acceleration storage.

use crate::{GRAVITY, particle::ParticleState, utils::VectorSeries};

/// Per-particle acceleration components used by the integrator.
///
/// The three vectors are parallel to the vectors in [`ParticleState`]: entry
/// `i` in each vector belongs to the particle at index `i`.
pub struct ForceBuffer {
    // /// X acceleration components in AU/year².
    // pub ax: Vec<f64>,
    // /// Y acceleration components in AU/year².
    // pub ay: Vec<f64>,
    // /// Z acceleration components in AU/year².
    // pub az: Vec<f64>,
    pub acceleration: VectorSeries,
}

impl ForceBuffer {
    /// Creates a zeroed acceleration buffer for `number_particles` particles.
    pub fn new(number_particles: usize) -> Self {
        ForceBuffer {
            acceleration: VectorSeries {
                x: vec![0.0; number_particles],
                y: vec![0.0; number_particles],
                z: vec![0.0; number_particles],
            },
        }
    }

    /// Recomputes and stores the acceleration of every particle.
    ///
    /// Massive particles contribute to the acceleration of every other
    /// particle. Massless particles still receive acceleration but do not
    /// contribute to the force calculation. Self-interaction is skipped.
    /// Existing values in this buffer are overwritten.
    pub fn update_accelerations(&mut self, state: &ParticleState) {
        for target_index in 0..state.mass.len() {
            if !state.alive[target_index] {
                continue;
            }

            let mut ax = 0.0;
            let mut ay = 0.0;
            let mut az = 0.0;

            for source_index in 0..state.mass.len() {
                if target_index == source_index {
                    continue;
                }

                if !state.alive[source_index] {
                    continue;
                }

                let Some(source_mass) = state.mass[source_index] else {
                    continue;
                };

                let dx = state.position.x[target_index] - state.position.x[source_index];
                let dy = state.position.y[target_index] - state.position.y[source_index];
                let dz = state.position.z[target_index] - state.position.z[source_index];

                let dist_squared = dx * dx + dy * dy + dz * dz;
                ax += gravity_acceleration(dx, dist_squared, source_mass);
                ay += gravity_acceleration(dy, dist_squared, source_mass);
                az += gravity_acceleration(dz, dist_squared, source_mass);
            }

            self.acceleration.x[target_index] = ax;
            self.acceleration.y[target_index] = ay;
            self.acceleration.z[target_index] = az;
        }
    }
}

/// Computes one Cartesian component of gravitational acceleration.
///
/// `dimension_dist` is the target coordinate minus the source coordinate,
/// `dist_squared` is the full three-dimensional separation squared, and
/// `attractor_mass` is the source mass in solar-mass units.
///
/// The calculation is singular when `dist_squared` is zero; callers should
/// prevent coincident source and target positions when appropriate.
pub fn gravity_acceleration(dimension_dist: f64, dist_squared: f64, attractor_mass: f64) -> f64 {
    // fvec = m1 avec = g m1 m2 / rmag^3 rvec
    -GRAVITY * attractor_mass * dimension_dist / (dist_squared * dist_squared.sqrt())
}
