//! `nsim` is a small Rust N-body simulation prototype. It currently models
//! massive bodies and massless test particles using a fixed-timestep
//! leapfrog/velocity-Verlet integrator.
#![warn(missing_docs)]
#![warn(clippy::all, clippy::pedantic)]
#![warn(rustdoc::broken_intra_doc_links)]

use std::f64::consts::PI;

use nsim::{
    diagnostics::Diagnostics,
    force::{ForceSystem, GRAVITY, NewtonianGravity},
    integration::leapfrog_timestep,
    math_util::vector3::Vector3,
    output::{append_particle_timestep, create_particle_file, write_diagnostics_file},
    particle::{Particle, ParticleSystem},
};

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
        position: Vector3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        velocity: Vector3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        mass: 1.0,
    });

    system.new_particle(Particle {
        name: String::from("Jupiter"),
        radius: 1.0,
        position: Vector3 {
            x: 5.0,
            y: 0.0,
            z: 0.0,
        },
        velocity: Vector3 {
            x: 0.0,
            y: f64::sqrt(GRAVITY / 5.0),
            z: 0.0,
        },
        mass: 0.0,
    });

    // ------------------------------------------------------------------------

    let mut time = time_start;
    let mut forces = ForceSystem::new(system.particle_count());
    forces.add_pairwise_force(NewtonianGravity);
    let mut diagnostics = Diagnostics::default();

    // --------------------- RECORD INITIAL STATE -----------------------------
    for particle_index in 0..system.particle_count() {
        create_particle_file(&system, particle_index)?;
        append_particle_timestep(&system, particle_index, time_start)?;
    }

    let initial_evaluation = forces.evaluate(system.state());
    diagnostics.record(
        time,
        system.state(),
        initial_evaluation.potential_energy.total(),
    );
    // ------------------------------------------------------------------------

    while time <= time_end {
        let force_evaluation = leapfrog_timestep(system.state_mut(), &mut forces, time_step);

        time += time_step;

        for particle_index in 0..system.particle_count() {
            append_particle_timestep(&system, particle_index, time)?;
        }

        diagnostics.record(
            time,
            system.state(),
            force_evaluation.potential_energy.total(),
        );
    }

    write_diagnostics_file(&diagnostics)?;

    Ok(())
}
