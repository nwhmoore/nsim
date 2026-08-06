//! `nsim` is a small Rust N-body simulation prototype. It currently models
//! massive bodies and massless test particles using a fixed-timestep
//! leapfrog/velocity-Verlet integrator.
#![warn(missing_docs)]
#![warn(clippy::all, clippy::pedantic)]
#![warn(rustdoc::broken_intra_doc_links)]

use std::f64::consts::PI;

use crate::{
    diagnostics::Diagnostics,
    force::ForceBuffer,
    integration::leapfrog_timestep,
    output::{append_particle_timestep, create_particle_file, write_diagnostics_file},
    particle::{Particle, ParticleSystem},
};

mod diagnostics;
mod force;
mod integration;
mod output;
mod particle;
mod utils;

/// Gravitational constant in AU³ · year⁻² · solar-mass⁻¹.
///
/// The units of this constant set the units of the entire simulation.
const GRAVITY: f64 = 4.0 * PI * PI;

fn main() -> std::io::Result<()> {
    // ---------------------------  INITIAL PARAMETERS ------------------------

    let time_start = 0.0;
    let time_end = 2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt();
    // 1% of period
    let time_step = time_end * 0.01;

    let mut system = ParticleSystem::new_system();

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

    // ------------------------------------------------------------------------

    let mut time = time_start;
    let mut force_buffer = ForceBuffer::new(system.catalog.id.len());
    let mut diagnostics = Diagnostics::default();

    // --------------------- RECORD INITIAL STATE -----------------------------
    for particle_index in 0..system.catalog.name.len() {
        create_particle_file(&system, particle_index)?;
        append_particle_timestep(&system, particle_index, time_start)?;
    }

    let initial_evaluation = force_buffer.compute_accelerations(&system.state);
    diagnostics.record(time, &system.state, initial_evaluation.potential_energy);
    // ------------------------------------------------------------------------

    while time <= time_end {
        let force_evaluation = leapfrog_timestep(&mut system.state, &mut force_buffer, time_step);

        time += time_step;

        for particle_index in 0..system.catalog.name.len() {
            append_particle_timestep(&system, particle_index, time)?;
        }

        diagnostics.record(time, &system.state, force_evaluation.potential_energy);
    }

    write_diagnostics_file(&diagnostics)?;

    Ok(())
}
