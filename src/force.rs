//! Acceleration calculation and reusable storage.

use crate::{
    force::gravity::{
        GRAVITY, gravitational_potential_energy, gravity_acceleration,
    }, particle::ParticleState, utils::{Geometry, KahanAccumulator, Vector3, Vector3Series},
};

pub mod gravity;

/// Quantities calculated alongside one gravitational force evaluation.
pub struct ForceEvaluation {
    /// Gravitational potential energy of the active massive bodies.
    pub potential_energy: f64,
}

/// Per-particle acceleration components used by the integrator.
///
/// The three vectors are parallel to the vectors in [`ParticleState`]: entry
/// `i` in each vector belongs to the particle at index `i`.
pub struct ForceBuffer {
    /// Cartesian acceleration components.
    ///
    /// The component vectors are aligned with the particle indices in the
    /// [`ParticleState`] used for the most recent evaluation.
    pub acceleration: Vector3Series,
    accumulator: AccelerationAccumulator,
    active_massive: Vec<(usize, f64)>,
    active_massless: Vec<usize>,
}

/// Per-particle Kahan-compensated acceleration totals used during one force
/// evaluation.
struct AccelerationAccumulator {
    /// Accumulated X-component accelerations for each particle.
    x: Vec<KahanAccumulator>,
    /// Accumulated Y-component accelerations for each particle.
    y: Vec<KahanAccumulator>,
    /// Accumulated Z-component accelerations for each particle.
    z: Vec<KahanAccumulator>,
}

impl AccelerationAccumulator {
    /// Creates an accumulator with one compensated total per particle.
    fn new(number_particles: usize) -> Self {
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
    fn add(&mut self, particle_idx: usize, acceleration: &Vector3) {
        self.x[particle_idx].add(acceleration.x);
        self.y[particle_idx].add(acceleration.y);
        self.z[particle_idx].add(acceleration.z);
    }
}

impl ForceBuffer {
    /// Creates a zeroed acceleration buffer for `number_particles` particles.
    pub fn new(number_particles: usize) -> Self {
        ForceBuffer {
            acceleration: Vector3Series {
                x: vec![0.0; number_particles],
                y: vec![0.0; number_particles],
                z: vec![0.0; number_particles],
            },
            accumulator: AccelerationAccumulator::new(number_particles),
            active_massive: Vec::with_capacity(number_particles),
            active_massless: Vec::with_capacity(number_particles),
        }
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
        // Reuse the preallocated classification buffers for this evaluation.
        self.fill_alive_particle_groups(state);

        for &particle_idx in self
            .active_massive
            .iter()
            .map(|(particle_idx, _)| particle_idx)
            .chain(self.active_massless.iter())
        {
            self.accumulator.x[particle_idx].reset();
            self.accumulator.y[particle_idx].reset();
            self.accumulator.z[particle_idx].reset();
        }

        // Pairwise forces
        let mut grav_pot_ener = KahanAccumulator::default();

        for (i, &(first_idx, first_mass)) in self.active_massive.iter().enumerate() {
            for &(second_idx, second_mass) in &self.active_massive[i + 1..] {
                let geometry = Geometry::calculate_geometry(Vector3 {
                    x: state.position.x[first_idx] - state.position.x[second_idx],
                    y: state.position.y[first_idx] - state.position.y[second_idx],
                    z: state.position.z[first_idx] - state.position.z[second_idx],
                });

                // Gravity
                let scale = -GRAVITY * geometry.inv_dist_cubed;
                let first_acceleration = gravity_acceleration(second_mass, &geometry.r_vec, scale);
                let second_acceleration = gravity_acceleration(first_mass, &geometry.r_vec, scale);
                self.accumulator.add(first_idx, &first_acceleration);
                self.accumulator.add(second_idx, &second_acceleration);

                grav_pot_ener.add(gravitational_potential_energy(
                    first_mass,
                    second_mass,
                    &geometry,
                ));
            }

            self.acceleration.x[first_idx] = self.accumulator.x[first_idx].total();
            self.acceleration.y[first_idx] = self.accumulator.y[first_idx].total();
            self.acceleration.z[first_idx] = self.accumulator.z[first_idx].total();
        }

        // One-way interactions for massless particles
        for &small_particle_idx in &self.active_massless {
            for &(large_particle_idx, attractor_mass) in &self.active_massive {
                // Geometry
                let geometry = Geometry::calculate_geometry(Vector3 {
                    x: state.position.x[small_particle_idx] - state.position.x[large_particle_idx],
                    y: state.position.y[small_particle_idx] - state.position.y[large_particle_idx],
                    z: state.position.z[small_particle_idx] - state.position.z[large_particle_idx],
                });

                // Gravity
                let scale = -GRAVITY * geometry.inv_dist_cubed;
                let gravity_acceleration = gravity_acceleration(attractor_mass, &geometry.r_vec, scale);
                self.accumulator
                    .add(small_particle_idx, &gravity_acceleration);
            }
            self.acceleration.x[small_particle_idx] =
                self.accumulator.x[small_particle_idx].total();
            self.acceleration.y[small_particle_idx] =
                self.accumulator.y[small_particle_idx].total();
            self.acceleration.z[small_particle_idx] =
                self.accumulator.z[small_particle_idx].total();
        }

        ForceEvaluation {
            potential_energy: grav_pot_ener.total(),
        }
    }

    /// Reclassifies active particles into massive and massless groups for the
    /// current force evaluation.
    fn fill_alive_particle_groups(&mut self, state: &ParticleState) {
        self.active_massive.clear();
        self.active_massless.clear();

        for particle_idx in 0..state.alive.len() {
            if !state.alive[particle_idx] {
                continue;
            }

            if let Some(mass) = state.mass[particle_idx] {
                self.active_massive.push((particle_idx, mass));
            } else {
                self.active_massless.push(particle_idx);
            }
        }
    }
}
