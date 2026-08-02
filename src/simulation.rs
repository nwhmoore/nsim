use crate::{
    body::{LargeBody, Positioned, SmallBody},
    math::{Acceleration, gravity_acceleration, time::Time},
};

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
