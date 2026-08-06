//! Fixed-timestep leapfrog/velocity-Verlet integration.
//!
//! The integrator reuses the acceleration already stored in [`ForceBuffer`]
//! for the current state. The caller must perform one initial force evaluation
//! before the first timestep. Each completed timestep leaves the buffer ready
//! for the next timestep.

use crate::{
    force::{ForceBuffer, ForceEvaluation},
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
    force_buffer: &mut ForceBuffer,
    fixed_dt: f64,
) -> ForceEvaluation {
    for particle_index in 0..state.mass.len() {
        if !state.alive[particle_index] {
            continue;
        }

        state.velocity.x[particle_index] +=
            force_buffer.acceleration.x[particle_index] * 0.5 * fixed_dt;
        state.position.x[particle_index] += state.velocity.x[particle_index] * fixed_dt;

        state.velocity.y[particle_index] +=
            force_buffer.acceleration.y[particle_index] * 0.5 * fixed_dt;
        state.position.y[particle_index] += state.velocity.y[particle_index] * fixed_dt;

        state.velocity.z[particle_index] +=
            force_buffer.acceleration.z[particle_index] * 0.5 * fixed_dt;
        state.position.z[particle_index] += state.velocity.z[particle_index] * fixed_dt;
    }

    let force_evaluation = force_buffer.compute_accelerations(state);

    for particle_index in 0..state.mass.len() {
        if !state.alive[particle_index] {
            continue;
        }

        state.velocity.x[particle_index] +=
            force_buffer.acceleration.x[particle_index] * 0.5 * fixed_dt;
        state.velocity.y[particle_index] +=
            force_buffer.acceleration.y[particle_index] * 0.5 * fixed_dt;
        state.velocity.z[particle_index] +=
            force_buffer.acceleration.z[particle_index] * 0.5 * fixed_dt;
    }

    force_evaluation
}
