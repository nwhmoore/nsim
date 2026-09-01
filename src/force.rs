//! Acceleration calculation and reusable storage.

use crate::{math_util::Vector3Series, particle::ParticleState};

mod constant;
mod drag;
mod gravity;
mod harmonic;

pub use constant::*;
pub use drag::*;
pub use gravity::*;
pub use harmonic::*;

/// forces being added to the simulation
#[derive(Default, Clone)]
pub struct ForceConfiguration {
    forces: Vec<Box<dyn Force>>,
}

impl ForceConfiguration {
    /// add a force to the simulation
    ///
    /// note: do not add a single force more than once
    pub fn add_force<F>(&mut self, force: F)
    where
        F: Force + 'static,
    {
        self.forces.push(Box::new(force));
    }
}

/// force system which holds the config and the internal acceleration buffer for
/// all particles
pub struct ForceSystem {
    configuration: ForceConfiguration,
    buffer: ForceBuffer,
}

impl ForceSystem {
    /// creates a new force system
    #[must_use]
    pub fn new(configuration: ForceConfiguration, particle_count: usize) -> Self {
        Self {
            configuration,
            buffer: ForceBuffer::new(particle_count),
        }
    }

    /// evaluates all forces on all particles in a given state
    pub fn evaluate(&mut self, particle_state: &ParticleState) {
        self.buffer.clear();

        let mut output = ForceEvaluation {
            accelerations: &mut self.buffer.accelerations,
        };

        for force in &self.configuration.forces {
            force.evaluate(particle_state, &mut output);
        }
    }

    /// returns the force buffer
    #[must_use]
    pub fn buffer(&self) -> &ForceBuffer {
        &self.buffer
    }

    /// returns the configured forces
    #[must_use]
    pub fn configured_forces(&self) -> &[Box<dyn Force>] {
        &self.configuration.forces
    }
}

/// defines a force
pub trait Force: ForceClone {
    /// evaluates this force on all particles
    fn evaluate(&self, particle_state: &ParticleState, output: &mut ForceEvaluation<'_>);

    /// calculates the potential energy associated with this force
    fn calculate_potential_energy(&self, _particle_state: &ParticleState) -> Option<f64> {
        None
    }
}

/// scratch work space for the force evaluation
pub struct ForceEvaluation<'a> {
    accelerations: &'a mut Vector3Series,
}

/// buffered accelerations of the particles
pub struct ForceBuffer {
    accelerations: Vector3Series,
}

impl ForceBuffer {
    /// Creates a zeroed acceleration buffer for `number_particles` particles.
    fn new(number_particles: usize) -> Self {
        ForceBuffer {
            accelerations: Vector3Series::new_with_zeros(number_particles),
        }
    }

    fn clear(&mut self) {
        self.accelerations.x.fill(0.0);
        self.accelerations.y.fill(0.0);
        self.accelerations.z.fill(0.0);
    }

    /// Returns the per-particle acceleration vectors stored by the buffer.
    #[must_use]
    pub fn accelerations(&self) -> &Vector3Series {
        &self.accelerations
    }
}

/// lets us clone the force system
pub trait ForceClone {
    /// clones
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
