//! `nsim` is a small Rust N-body simulation prototype. It currently models
//! massive bodies and massless test particles using a fixed-timestep
//! leapfrog/velocity-Verlet integrator.
#![warn(missing_docs)]
#![warn(clippy::all, clippy::pedantic)]
#![warn(rustdoc::broken_intra_doc_links)]

use std::f64::consts::PI;

use nsim::{
    force::{GRAVITY, NewtonianGravity},
    integration::Leapfrog,
    math_util::Vector3,
    particle::{Particle, ParticleSystem},
    simulation::Simulation,
};

fn main() {
    // ---------------------------  INITIAL PARAMETERS ------------------------
    let mut particle_system = ParticleSystem::new_system();

    particle_system.add_particle(Particle {
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

    particle_system.add_particle(Particle {
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

    // one period
    let end_time = 2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt();

    let mut simulation = Simulation::new()
        .with_particle_system(particle_system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .set_time_step(end_time * 0.01)
        .set_diagnostic_interval(end_time)
        .build();

    simulation.run_until(end_time);
}
