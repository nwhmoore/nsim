//! Fixed-timestep leapfrog/velocity-Verlet integration.
//!
//! The integrator reuses the acceleration already stored in [`ForceBuffer`]
//! for the current state. The caller must perform one initial force evaluation
//! before the first timestep. Each completed timestep leaves the buffer ready
//! for the next timestep.

use crate::{force::ForceSystem, particle::ParticleState};

mod leapfrog;
mod runge_kutta;

pub use leapfrog::*;
pub use runge_kutta::*;

/// which integrator that advances the simulation
pub trait Integrator {
    /// performs any required initialization when simulation building
    fn initialize(&mut self, particle_state: &ParticleState);

    /// advances the simulation one timestep
    fn evaluate_timestep(&mut self, state: &mut ParticleState, forces: &mut ForceSystem, dt: f64);

    /// provides a warning when simulation building
    fn warn();
}

/// empy integrator used as a default or used for testing.
pub struct NoIntegrator;

impl Integrator for NoIntegrator {
    fn initialize(&mut self, _state: &ParticleState) {}

    fn evaluate_timestep(
        &mut self,
        _state: &mut ParticleState,
        _forces: &mut ForceSystem,
        _dt: f64,
    ) {
    }

    fn warn() {
        eprintln!("[WARNING] :: no integrator selected; simulation will not advance.");
    }
}
