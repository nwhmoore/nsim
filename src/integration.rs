//! Fixed-timestep leapfrog/velocity-Verlet integration.
//!
//! The integrator reuses the acceleration already stored in [`ForceBuffer`]
//! for the current state. The caller must perform one initial force evaluation
//! before the first timestep. Each completed timestep leaves the buffer ready
//! for the next timestep.

use crate::{
    force::{ForceDiagnostics, ForceSystem},
    particle::ParticleState,
};

pub trait Integrator {
    fn evaluate_timestep(
        &mut self,
        state: &mut ParticleState,
        forces: &mut ForceSystem,
        dt: f64,
    ) -> ForceDiagnostics;
}

#[derive(Clone)]
pub struct Leapfrog;

impl Integrator for Leapfrog {
    /// Advances all bodies by one timestep using kick-drift-kick integration.
    ///
    /// The update sequence is:
    ///
    /// 1. Use the cached acceleration at the current positions for the first kick.
    /// 2. Advance velocities by half a timestep.
    /// 3. Advance positions by one full timestep using the half-step velocities.
    /// 4. Recompute acceleration and potential energy at the new positions.
    /// 5. Advance velocities by the remaining half timestep.
    ///
    /// The state and force-buffer vectors must have matching lengths. Before the
    /// call, the force buffer must contain accelerations evaluated at the input
    /// positions. After the call, it contains accelerations evaluated at the
    /// output positions, and can be reused by the next call.
    ///
    /// # Returns
    ///
    /// The force evaluation at the output state, including its gravitational
    /// potential energy.
    ///
    /// # Panics
    ///
    /// This function may panic if the state and force-buffer component vectors do
    /// not have matching lengths.
    fn evaluate_timestep(
        &mut self,
        state: &mut ParticleState,
        forces: &mut ForceSystem,
        dt: f64,
    ) -> ForceDiagnostics {
        let n = state.particle_count();
        let half_dt = 0.5 * dt;

        {
            let accelerations = forces.buffer().accelerations();
            let (positions, velocities) = state.positions_and_velocities_mut();

            for i in 0..n {
                velocities.x[i] += accelerations.x[i] * half_dt;
                velocities.y[i] += accelerations.y[i] * half_dt;
                velocities.z[i] += accelerations.z[i] * half_dt;

                positions.x[i] += velocities.x[i] * dt;
                positions.y[i] += velocities.y[i] * dt;
                positions.z[i] += velocities.z[i] * dt;
            }
        }

        let force_evaluation = forces.evaluate(state);

        {
            let accelerations = forces.buffer().accelerations();
            let velocities = state.velocities_mut();

            for i in 0..n {
                velocities.x[i] += accelerations.x[i] * half_dt;
                velocities.y[i] += accelerations.y[i] * half_dt;
                velocities.z[i] += accelerations.z[i] * half_dt;
            }
        }

        force_evaluation
    }
}

pub struct NoIntegrator;

impl Integrator for NoIntegrator {
    fn evaluate_timestep(
        &mut self,
        _state: &mut ParticleState,
        _forces: &mut ForceSystem,
        _dt: f64,
    ) -> ForceDiagnostics {
        ForceDiagnostics::default()
    }
}
