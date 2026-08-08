//! Gravitational acceleration, potential-energy calculation, and reusable
//! acceleration storage.

use crate::{
    GRAVITY,
    particle::ParticleState,
    utils::{KahanAccumulator, VectorSeries},
};

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
    pub acceleration: VectorSeries,
    accumulator: AccelerationAccumulator,
    active_massive: Vec<(usize, f64)>,
    active_massless: Vec<usize>,
}

struct AccelerationAccumulator {
    x: Vec<KahanAccumulator>,
    y: Vec<KahanAccumulator>,
    z: Vec<KahanAccumulator>,
}

impl AccelerationAccumulator {
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
        // todo: preallocate these vecs
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
                let dx = state.position.x[first_idx] - state.position.x[second_idx];
                let dy = state.position.y[first_idx] - state.position.y[second_idx];
                let dz = state.position.z[first_idx] - state.position.z[second_idx];

                let dist_squared = dx * dx + dy * dy + dz * dz;
                let dist = dist_squared.sqrt();
                let inv_dist_cubed = 1.0 / (dist * dist_squared);
                let acceleration_scale_from_first = -GRAVITY * first_mass * inv_dist_cubed;
                let acceleration_scale_from_second = -GRAVITY * second_mass * inv_dist_cubed;

                self.accumulator.x[first_idx].add(dx * acceleration_scale_from_second);
                self.accumulator.y[first_idx].add(dy * acceleration_scale_from_second);
                self.accumulator.z[first_idx].add(dz * acceleration_scale_from_second);

                self.accumulator.x[second_idx].add(-dx * acceleration_scale_from_first);
                self.accumulator.y[second_idx].add(-dy * acceleration_scale_from_first);
                self.accumulator.z[second_idx].add(-dz * acceleration_scale_from_first);

                grav_pot_ener.add(gravity_potential(first_mass, second_mass, dist));
            }

            self.acceleration.x[first_idx] = self.accumulator.x[first_idx].total();
            self.acceleration.y[first_idx] = self.accumulator.y[first_idx].total();
            self.acceleration.z[first_idx] = self.accumulator.z[first_idx].total();
        }

        // massless particle interactions
        // TODO: apply the same acceleration scale from attractor approach here.
        for &small_particle_idx in &self.active_massless {
            for &(large_particle_idx, attractor_mass) in &self.active_massive {
                let dx =
                    state.position.x[small_particle_idx] - state.position.x[large_particle_idx];
                let dy =
                    state.position.y[small_particle_idx] - state.position.y[large_particle_idx];
                let dz =
                    state.position.z[small_particle_idx] - state.position.z[large_particle_idx];

                let dist_squared = dx * dx + dy * dy + dz * dz;
                let dist = dist_squared.sqrt();
                let inv_dist_cubed = 1.0 / (dist * dist_squared);
                let acceleration_scale_from_attractor = -GRAVITY * attractor_mass * inv_dist_cubed;

                self.accumulator.x[small_particle_idx].add(dx * acceleration_scale_from_attractor);
                self.accumulator.y[small_particle_idx].add(dy * acceleration_scale_from_attractor);
                self.accumulator.z[small_particle_idx].add(dz * acceleration_scale_from_attractor);
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

/// Computes the Newtonian gravitational potential energy of one massive pair.
///
/// `dist_squared` is the squared separation of the bodies. The caller is
/// responsible for supplying nonzero separation and for counting each pair
/// only once.
#[must_use]
pub fn gravity_potential(first_mass: f64, second_mass: f64, dist: f64) -> f64 {
    -GRAVITY * first_mass * second_mass / dist
}
