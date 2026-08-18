//! Acceleration calculation and reusable storage.

use crate::{
    math_util::{
        Geometry,
        kahan::{Kahan3Series, KahanAccumulator},
        vector3::{Vector3, Vector3Series},
    },
    particle::ParticleState,
};

mod gravity;

pub use gravity::GRAVITY;
pub use gravity::NewtonianGravity;

pub struct ForceSystem<A = DirectAccumulator> {
    pairwise: Vec<Box<dyn PairwiseForce>>,
    accumulator: A,
    buffer: ForceBuffer,
}

impl<A: AccelerationAccumulator> ForceSystem<A> {
    pub fn new(particle_count: usize) -> Self {
        ForceSystem {
            pairwise: Vec::new(),
            accumulator: A::with_len(particle_count),
            buffer: ForceBuffer::new(particle_count),
        }
    }

    pub fn add_pairwise_force<F>(&mut self, force: F)
    where
        F: PairwiseForce + 'static,
    {
        self.pairwise.push(Box::new(force));
    }

    pub fn evaluate(&mut self, state: &ParticleState) -> ForceDiagnostics {
        let mut diagnostics = ForceDiagnostics::zero();

        self.accumulator.clear();

        for first_idx in 0..state.particle_count() {
            if !state.alive_statuses()[first_idx] {
                continue;
            }

            for second_idx in (first_idx + 1)..state.particle_count() {
                let geometry = Geometry::calculate_geometry(
                    state.positions().value_at(first_idx) - state.positions().value_at(second_idx),
                );

                for force in self.pairwise.iter() {
                    let contribution = force.evaluate_pair(state, first_idx, second_idx, &geometry);

                    self.accumulator
                        .add(first_idx, &contribution.first_acceleration);
                    self.accumulator
                        .add(second_idx, &contribution.second_acceleration);

                    diagnostics
                        .potential_energy
                        .add(contribution.potential_energy);
                }
            }
        }

        self.accumulator.finish(&mut self.buffer);

        diagnostics
    }

    pub fn buffer(&self) -> &ForceBuffer {
        &self.buffer
    }
}

pub trait PairwiseForce {
    fn evaluate_pair(
        &self,
        state: &ParticleState,
        first_idx: usize,
        second_idx: usize,
        geometry: &Geometry,
    ) -> PairForceContribution;
}

// In future, non-conservative forces may not apply to every force
pub struct PairForceContribution {
    pub first_acceleration: Vector3,
    pub second_acceleration: Vector3,
    pub potential_energy: f64,
}

pub trait AccelerationAccumulator {
    fn with_len(number_particles: usize) -> Self;

    fn clear(&mut self);

    fn add(&mut self, particle_idx: usize, acceleration: &Vector3);

    fn finish(&mut self, buffer: &mut ForceBuffer);
}

pub struct DirectAccumulator {
    accumulator: Vector3Series,
}

impl AccelerationAccumulator for DirectAccumulator {
    fn with_len(number_particles: usize) -> Self {
        DirectAccumulator {
            accumulator: Vector3Series::new_zeros(number_particles),
        }
    }

    fn clear(&mut self) {
        for idx in 0..self.accumulator.len() {
            self.accumulator.x[idx] = 0.0;
            self.accumulator.y[idx] = 0.0;
            self.accumulator.z[idx] = 0.0;
        }
    }

    fn add(&mut self, particle_idx: usize, acceleration: &Vector3) {
        self.accumulator.x[particle_idx] += acceleration.x;
        self.accumulator.y[particle_idx] += acceleration.y;
        self.accumulator.z[particle_idx] += acceleration.z;
    }

    fn finish(&mut self, buffer: &mut ForceBuffer) {
        for idx in 0..self.accumulator.len() {
            buffer.accelerations.x[idx] = self.accumulator.x[idx];
            buffer.accelerations.y[idx] = self.accumulator.y[idx];
            buffer.accelerations.z[idx] = self.accumulator.z[idx];
        }
    }
}

pub struct CompensatedAccumulator {
    accumulator: Kahan3Series,
}

impl AccelerationAccumulator for CompensatedAccumulator {
    fn with_len(number_particles: usize) -> Self {
        CompensatedAccumulator {
            accumulator: Kahan3Series::new(number_particles),
        }
    }

    fn clear(&mut self) {
        for idx in 0..self.accumulator.len() {
            self.accumulator.reset_at(idx);
        }
    }

    fn add(&mut self, particle_idx: usize, acceleration: &Vector3) {
        self.accumulator.add(particle_idx, acceleration)
    }

    fn finish(&mut self, buffer: &mut ForceBuffer) {
        for idx in 0..self.accumulator.len() {
            buffer
                .accelerations
                .set_value_at(idx, self.accumulator.total(idx));
        }
    }
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
    accelerations: Vector3Series,
}

impl ForceBuffer {
    /// Creates a zeroed acceleration buffer for `number_particles` particles.
    fn new(number_particles: usize) -> Self {
        ForceBuffer {
            accelerations: Vector3Series::new_zeros(number_particles),
        }
    }

    /// Returns the per-particle acceleration vectors stored by the buffer.
    pub fn accelerations(&self) -> &Vector3Series {
        &self.accelerations
    }
}

/// Quantities calculated alongside one force evaluation.
pub struct ForceDiagnostics {
    /// Gravitational potential energy of the active massive bodies.
    pub potential_energy: KahanAccumulator,
}

impl ForceDiagnostics {
    pub fn zero() -> Self {
        ForceDiagnostics {
            potential_energy: KahanAccumulator::default(),
        }
    }
}
