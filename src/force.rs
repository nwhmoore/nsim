//! Gravitational acceleration calculation and reusable acceleration storage.

use crate::{
    GRAVITY,
    particle::ParticleState,
    utils::{KahanAccumulator, VectorSeries},
};

pub struct ForceEvaluation {
    pub potential_energy: f64,
}

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
    pub fn compute_accelerations(&mut self, state: &ParticleState) -> ForceEvaluation {
        let mut grav_pot_ener = KahanAccumulator::default();
        for target_index in 0..state.mass.len() {
            if !state.alive[target_index] {
                continue;
            }

            let target_mass = state.mass[target_index];

            let mut ax = KahanAccumulator::default();
            let mut ay = KahanAccumulator::default();
            let mut az = KahanAccumulator::default();

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

                ax.add(gravity_acceleration(source_mass, dx, dist_squared));
                ay.add(gravity_acceleration(source_mass, dy, dist_squared));
                az.add(gravity_acceleration(source_mass, dz, dist_squared));

                if target_index < source_index
                    && let Some(target_mass) = target_mass
                {
                    grav_pot_ener.add(gravity_potential(target_mass, source_mass, dist_squared));
                }
            }

            self.acceleration.x[target_index] = ax.total();
            self.acceleration.y[target_index] = ay.total();
            self.acceleration.z[target_index] = az.total();
        }

        ForceEvaluation {
            potential_energy: grav_pot_ener.total(),
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
pub fn gravity_acceleration(attractor_mass: f64, dimension_dist: f64, dist_squared: f64) -> f64 {
    // fvec = m1 avec = g m1 m2 / rmag^3 rvec
    -GRAVITY * attractor_mass * dimension_dist / (dist_squared * dist_squared.sqrt())
}

pub fn gravity_potential(self_mass: f64, attractor_mass: f64, dist_squared: f64) -> f64 {
    -GRAVITY * self_mass * attractor_mass / dist_squared.sqrt()
}
