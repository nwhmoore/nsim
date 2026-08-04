//! Fixed-timestep leapfrog/velocity-Verlet integration.

use crate::{force::ForceBuffer, particle::ParticleState};

/// Advances all bodies by one timestep using kick-drift-kick integration.
///
/// The update sequence is:
///
/// 1. Compute accelerations at the current positions.
/// 2. Advance velocities by half a timestep.
/// 3. Advance positions by one full timestep using the half-step velocities.
/// 4. Recompute accelerations at the new positions.
/// 5. Advance velocities by the remaining half timestep.
pub fn leapfrog_timestep(state: &mut ParticleState, force_buffer: &mut ForceBuffer, dt: f64) {
    force_buffer.update_accelerations(state);

    // half-step velocity v_(t+.5dt) = v_t + .5 * a_t * dt

    for particle_index in 0..state.mass.len() {
        state.vx[particle_index] += force_buffer.ax[particle_index] * 0.5 * dt;
        state.x[particle_index] += state.vx[particle_index] * dt;

        state.vy[particle_index] += force_buffer.ay[particle_index] * 0.5 * dt;
        state.y[particle_index] += state.vy[particle_index] * dt;

        state.vz[particle_index] += force_buffer.az[particle_index] * 0.5 * dt;
        state.z[particle_index] += state.vz[particle_index] * dt;
    }

    // recompute accelerations at current positions
    force_buffer.update_accelerations(state);
    
    // finish velocity update
    for particle_index in 0..state.mass.len() {
        state.vx[particle_index] += force_buffer.ax[particle_index] * 0.5 * dt;
        state.vy[particle_index] += force_buffer.ay[particle_index] * 0.5 * dt;
        state.vz[particle_index] += force_buffer.az[particle_index] * 0.5 * dt;
    }
}
