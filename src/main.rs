use std::f64::consts::PI;

use crate::{
    body::{LargeBody, SmallBody},
    math::{Acceleration, GRAVITY, Position, Velocity, time::Time},
    output::{append_body_timestep, create_body_file},
    simulation::leapfrog_timestep,
};

mod body;
mod math;
mod output;
mod simulation;

fn main() -> std::io::Result<()> {
    // ---------------  INITIAL PARAMETERS -------------------

    // years
    let time_start = Time::years(0.0);
    // years
    let time_end = Time::years(2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt());
    // 1% of period
    let time_step = time_end * 0.01;

    let mut large_bodies = [LargeBody {
        name: String::from("Sol"),
        pos: Position {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        vel: Velocity {
            u: 0.0,
            v: 0.0,
            w: 0.0,
        },
        acc: Acceleration::default(),
        mass: 1.0,
    }];

    let mut small_bodies = [SmallBody {
        name: String::from("Jupiter"),
        pos: Position {
            x: 5.0,
            y: 0.0,
            z: 0.0,
        },
        vel: Velocity {
            u: 0.0,
            v: f64::sqrt(GRAVITY / 5.0),
            w: 0.0,
        },
        acc: Acceleration::default(),
    }];

    // -----------------------------------------------------------

    for body in &large_bodies {
        create_body_file(body)?;
        append_body_timestep(body, time_start)?;
    }

    for body in &small_bodies {
        create_body_file(body)?;
        append_body_timestep(body, time_start)?;
    }

    // integrate based on leapfrog / velocity verlet
    let mut time = time_start;
    let mut large_acceleration_scratch = vec![Acceleration::default(); large_bodies.len()];

    while time <= time_end {
        leapfrog_timestep(
            &mut small_bodies,
            &mut large_bodies,
            &mut large_acceleration_scratch,
            time_step,
        );

        time += time_step;

        for body in &large_bodies {
            append_body_timestep(body, time)?;
        }

        for body in &small_bodies {
            append_body_timestep(body, time)?;
        }
    }

    Ok(())
}
