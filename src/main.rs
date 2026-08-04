use std::f64::consts::PI;

use crate::{
    force::ForceBuffer,
    integration::leapfrog_timestep,
    output::{append_particle_timestep, create_particle_file},
    particle::{Particle, ParticleSystem},
};

mod force;
mod integration;
mod output;
mod particle;

// currently in units of AU^3 yr^-2 M_sol^-1
const GRAVITY: f64 = 4.0 * PI * PI;

fn main() -> std::io::Result<()> {
    // ---------------  INITIAL PARAMETERS -------------------

    // years
    let time_start = 0.0;
    // years
    let time_end = 2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt();
    // 1% of period
    let time_step = time_end * 0.01;

    let mut system = ParticleSystem::new();

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

    for particle_index in 0..system.catalog.name.len() {
        create_particle_file(&system, particle_index)?;
        append_particle_timestep(&system, particle_index, time_start)?;
    }

    let mut time = time_start;
    let mut force_buffer = ForceBuffer::new(system.catalog.id.len());

    while time <= time_end {
        leapfrog_timestep(&mut system.state, &mut force_buffer, time_step);

        time += time_step;

        for particle_index in 0..system.catalog.name.len() {
            append_particle_timestep(&system, particle_index, time)?;
        }
    }

    Ok(())
}
