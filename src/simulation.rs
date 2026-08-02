//! Fixed-timestep leapfrog/velocity-Verlet integration.

use crate::{
    body::{LargeBody, Positioned, SmallBody},
    math::{Acceleration, gravity_acceleration, time::Time},
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
/// All body slices are updated in place. `large_acceleration_scratch` must
/// have the same length as `large_bodies`; it stores large-body accelerations
/// while the source positions are still being read.
pub fn leapfrog_timestep(
    small_bodies: &mut [SmallBody],
    large_bodies: &mut [LargeBody],
    large_acceleration_scratch: &mut [Acceleration],
    time_step: Time,
) {
    assert_eq!(large_bodies.len(), large_acceleration_scratch.len());

    let half_time_step = time_step * 0.5;

    update_accelerations(small_bodies, large_bodies, large_acceleration_scratch);

    // half-step velocity v_(t+.5) = v_t + .5 * a_t * dt

    for body in small_bodies.iter_mut() {
        body.vel += body.acc * half_time_step;
    }

    for body in large_bodies.iter_mut() {
        body.vel += body.acc * half_time_step;
    }

    // update position with half velocity

    for body in small_bodies.iter_mut() {
        body.pos += body.vel * time_step;
    }

    for body in large_bodies.iter_mut() {
        body.pos += body.vel * time_step;
    }

    // recompute accelerations at current positions
    update_accelerations(small_bodies, large_bodies, large_acceleration_scratch);

    // finish velocity update

    for body in small_bodies.iter_mut() {
        body.vel += body.acc * half_time_step;
    }

    for body in large_bodies.iter_mut() {
        body.vel += body.acc * half_time_step;
    }
}

/// Recomputes the acceleration of every body at its current position.
///
/// Small bodies receive acceleration from every large body. Large bodies
/// receive acceleration from every other large body, with self-interaction
/// excluded. The scratch buffer must have one entry per large body.
pub fn update_accelerations(
    small_bodies: &mut [SmallBody],
    large_bodies: &mut [LargeBody],
    large_acceleration_scratch: &mut [Acceleration],
) {
    // compute accelerations at current positions

    for small_body in small_bodies.iter_mut() {
        small_body.acc = update_this_acceleration(small_body, large_bodies);
    }

    for (index, (next_acceleration, large_body)) in large_acceleration_scratch
        .iter_mut()
        .zip(large_bodies.iter())
        .enumerate()
    {
        let mut total = Acceleration::default();

        for (source_index, source_body) in large_bodies.iter().enumerate() {
            if source_index != index {
                total += gravity_acceleration(large_body, source_body);
            }
        }

        *next_acceleration = total;
    }
    for (large_body, acceleration) in large_bodies
        .iter_mut()
        .zip(large_acceleration_scratch.iter().copied())
    {
        large_body.acc = acceleration;
    }
}

/// Sums the acceleration contributions from all large bodies for one target.
///
/// The target may be any type implementing [`Positioned`], including either
/// `SmallBody` or `LargeBody`.
pub fn update_this_acceleration<B>(body: &B, all_large: &[LargeBody]) -> Acceleration
where
    B: Positioned,
{
    let mut total = Acceleration {
        r: 0.0,
        s: 0.0,
        t: 0.0,
    };

    for large in all_large {
        total += gravity_acceleration(body, large);
    }

    total
}
