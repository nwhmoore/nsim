//! Acceleration calculation and reusable storage.

use crate::{
    math_util::{Vector3Series, kahan::KahanAccumulator},
    particle::ParticleState,
};

mod gravity;
mod simple;

pub use gravity::*;
pub use simple::*;

#[derive(Default, Clone)]
pub struct ForceConfiguration {
    pairwise: Vec<Box<dyn Force>>,
}

impl ForceConfiguration {
    pub fn add_force<F>(&mut self, force: F)
    where
        F: Force + 'static,
    {
        self.pairwise.push(Box::new(force));
    }
}

pub struct ForceSystem {
    configuration: ForceConfiguration,
    buffer: ForceBuffer,
}

impl ForceSystem {
    pub fn new(configuration: ForceConfiguration, particle_count: usize) -> Self {
        Self {
            configuration,
            buffer: ForceBuffer::new(particle_count),
        }
    }

    pub fn evaluate(&mut self, particle_state: &ParticleState) -> ForceDiagnostics {
        self.buffer.clear();

        let mut output = ForceEvaluation {
            accelerations: &mut self.buffer.accelerations,
            potential_energy: KahanAccumulator::default(),
        };

        for force in &self.configuration.pairwise {
            force.evaluate(particle_state, &mut output);
        }

        ForceDiagnostics {
            potential_energy: output.potential_energy.total(),
        }
    }

    pub fn buffer(&self) -> &ForceBuffer {
        &self.buffer
    }
}

pub trait Force: ForceClone {
    fn evaluate(&self, particle_state: &ParticleState, output: &mut ForceEvaluation<'_>);
}

pub struct ForceEvaluation<'a> {
    pub accelerations: &'a mut Vector3Series,
    pub potential_energy: KahanAccumulator,
}

pub struct ForceBuffer {
    pub accelerations: Vector3Series,
}

impl ForceBuffer {
    /// Creates a zeroed acceleration buffer for `number_particles` particles.
    fn new(number_particles: usize) -> Self {
        ForceBuffer {
            accelerations: Vector3Series::new_zeros(number_particles),
        }
    }

    fn clear(&mut self) {
        self.accelerations.x.fill(0.0);
        self.accelerations.y.fill(0.0);
        self.accelerations.z.fill(0.0);
    }

    /// Returns the per-particle acceleration vectors stored by the buffer.
    pub fn accelerations(&self) -> &Vector3Series {
        &self.accelerations
    }
}

/// Quantities calculated alongside one force evaluation.
#[derive(Default)]
pub struct ForceDiagnostics {
    /// Gravitational potential energy of the active massive bodies.
    pub potential_energy: f64,
}

pub trait ForceClone {
    fn clone_box(&self) -> Box<dyn Force>;
}

impl<T> ForceClone for T
where
    T: Force + Clone + 'static,
{
    fn clone_box(&self) -> Box<dyn Force> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Force> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}
