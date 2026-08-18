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
pub fn leapfrog_timestep(
    state: &mut ParticleState,
    forces: &mut ForceSystem,
    fixed_dt: f64,
) -> ForceDiagnostics {
    for particle_index in 0..state.particle_count() {
        if !state.alive_statuses()[particle_index] {
            continue;
        }

        let new_velocity = state.velocities().value_at(particle_index)
            + forces.buffer().accelerations().value_at(particle_index) * 0.5 * fixed_dt;
        state
            .velocities_mut()
            .set_value_at(particle_index, new_velocity);

        let new_position = state.positions().value_at(particle_index)
            + state.velocities().value_at(particle_index) * fixed_dt;
        state
            .positions_mut()
            .set_value_at(particle_index, new_position);
    }

    let force_evaluation = forces.evaluate(state);

    for particle_index in 0..state.particle_count() {
        if !state.alive_statuses()[particle_index] {
            continue;
        }

        let new_velocity = state.velocities().value_at(particle_index)
            + forces.buffer().accelerations().value_at(particle_index) * 0.5 * fixed_dt;
        state
            .velocities_mut()
            .set_value_at(particle_index, new_velocity);
    }

    force_evaluation
}
