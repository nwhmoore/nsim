//! Acceleration calculation and reusable storage.

use crate::{
    force::gravity::{GRAVITY, gravitational_potential_energy, gravity_acceleration},
    math_util::{
        Geometry,
        kahan::{Kahan3Series, KahanAccumulator},
        vector3::Vector3Series,
    },
    particle::ParticleState,
};

pub mod gravity;

/// Per-particle acceleration components used by the integrator.
///
/// The three vectors are parallel to the vectors in [`ParticleState`]: entry
/// `i` in each vector belongs to the particle at index `i`.
pub struct ForceBuffer {
    /// Cartesian acceleration components.
    ///
    /// The component vectors are aligned with the particle indices in the
    /// [`ParticleState`] used for the most recent evaluation.
    accelerations: Vector3Series,
    accumulator: Kahan3Series,
}

impl ForceBuffer {
    /// Creates a zeroed acceleration buffer for `number_particles` particles.
    pub fn new(number_particles: usize) -> Self {
        ForceBuffer {
            accelerations: Vector3Series::new_zeros(number_particles),
            accumulator: Kahan3Series::new(number_particles),
        }
    }

    /// Returns the per-particle acceleration vectors stored by the buffer.
    pub fn accelerations(&self) -> &Vector3Series {
        &self.accelerations
    }

    /// Recomputes and stores the acceleration of every active particle.
    ///
    /// Massive particles contribute to the acceleration of every other
    /// particle. Massless particles still receive acceleration but do not
    /// contribute to the force calculation. Self-interaction is skipped. The
    /// returned potential energy includes each pair of active massive bodies
    /// once. Existing acceleration values for active particles are overwritten;
    /// values for inactive particles are left unchanged.
    ///
    /// The position, velocity, mass, and activity vectors in `state` must have
    /// matching lengths. Coincident active particles produce the singular
    /// Newtonian force and potential and should be prevented by the caller.
    #[must_use]
    pub fn compute_accelerations(&mut self, state: &ParticleState) -> ForceEvaluation {
        for (particle_idx, _) in state.masses().iter().enumerate() {
            self.accumulator.reset_at(particle_idx);
        }

        // Pairwise forces
        let mut grav_pot_ener = KahanAccumulator::default();

        for (first_idx, &first_mass) in state.masses().iter().enumerate() {
            for (second_idx, &second_mass) in state.masses().iter().enumerate().skip(first_idx + 1)
            {
                let geometry = Geometry::calculate_geometry(
                    state.positions().value_at(first_idx) - state.positions().value_at(second_idx),
                );

                // Gravity
                let scale = -GRAVITY * geometry.inv_dist_cubed();
                let first_acceleration = gravity_acceleration(second_mass, geometry.r_vec(), scale);
                let second_acceleration =
                    gravity_acceleration(first_mass, geometry.r_vec(), -scale);
                self.accumulator.add(first_idx, &first_acceleration);
                self.accumulator.add(second_idx, &second_acceleration);

                grav_pot_ener.add(gravitational_potential_energy(
                    first_mass,
                    second_mass,
                    &geometry,
                ));
            }

            self.accelerations
                .set_value_at(first_idx, self.accumulator.total(first_idx));
        }

        ForceEvaluation {
            potential_energy: grav_pot_ener.total(),
        }
    }
}

/// Quantities calculated alongside one force evaluation.
pub struct ForceEvaluation {
    /// Gravitational potential energy of the active massive bodies.
    pub potential_energy: f64,
}
