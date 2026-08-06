//! Fixed-timestep leapfrog/velocity-Verlet integration.

use crate::{
    force::{ForceBuffer, ForceEvaluation},
    particle::ParticleState,
};

/// Advances all bodies by one timestep using kick-drift-kick integration.
///
/// The update sequence is:
///
/// 1. Compute accelerations at the current positions.
/// 2. Advance velocities by half a timestep.
/// 3. Advance positions by one full timestep using the half-step velocities.
/// 4. Recompute accelerations at the new positions.
/// 5. Advance velocities by the remaining half timestep.
///
/// The state and force-buffer vectors must have matching lengths and the force
/// buffer is reused in place.
pub fn leapfrog_timestep(
    state: &mut ParticleState,
    force_buffer: &mut ForceBuffer,
    fixed_dt: f64,
) -> ForceEvaluation {
    for particle_index in 0..state.mass.len() {
        if !state.alive[particle_index] {
            continue;
        }

        state.velocity.x[particle_index] += force_buffer.acceleration.x[particle_index] * 0.5 * fixed_dt;
        state.position.x[particle_index] += state.velocity.x[particle_index] * fixed_dt;

        state.velocity.y[particle_index] += force_buffer.acceleration.y[particle_index] * 0.5 * fixed_dt;
        state.position.y[particle_index] += state.velocity.y[particle_index] * fixed_dt;

        state.velocity.z[particle_index] += force_buffer.acceleration.z[particle_index] * 0.5 * fixed_dt;
        state.position.z[particle_index] += state.velocity.z[particle_index] * fixed_dt;
    }

    let force_evaluation = force_buffer.compute_accelerations(state);

    for particle_index in 0..state.mass.len() {
        if !state.alive[particle_index] {
            continue;
        }

        state.velocity.x[particle_index] += force_buffer.acceleration.x[particle_index] * 0.5 * fixed_dt;
        state.velocity.y[particle_index] += force_buffer.acceleration.y[particle_index] * 0.5 * fixed_dt;
        state.velocity.z[particle_index] += force_buffer.acceleration.z[particle_index] * 0.5 * fixed_dt;
    }

    force_evaluation
}
