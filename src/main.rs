use std::f64::consts::PI;

use crate::{
    integration::leapfrog_timestep,
    particle::{Particle, ParticleSystem},
};

mod force;
mod integration;
mod particle;

const GRAVITY: f64 = 4.0 * PI * PI;

fn main() -> std::io::Result<()> {
    // ---------------  INITIAL PARAMETERS -------------------

    // years
    let time_start = 0.0;
    // years
    let time_end = 2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt();
    // 1% of period
    let time_step = time_end * 0.01;

    let system = ParticleSystem::new();

    system.new_particle(Particle {
        name: String::from("Sol"),
        radius: 1.0,
        pos: (0.0, 0.0, 0.0),
        vel: (0.0, 0.0, 0.0),
        mass: Some(1.0),
    });

    system.new_particle(Particle {
        name: String::from("Jupiter"),
        radius: 1.0,
        pos: (5.0, 0.0, 0.0),
        vel: (0.0, f64::sqrt(GRAVITY / 5.0), 0.0),
        mass: None,
    });

    // -----------------------------------------------------------

    // integrate based on leapfrog / velocity verlet
    let mut time = time_start;

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
