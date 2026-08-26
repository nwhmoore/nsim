#![feature(test)]

use nsim::{
    force::NewtonianGravity,
    integration::Leapfrog,
    math_util::Vector3,
    particle::{Particle, ParticleSystem},
    simulation::SimulationBuilder,
};

fn solar_system(simulation_time: f64) {
    let mut particle_system = ParticleSystem::new_system();
    let velocity_scale = 365.2425;

    particle_system.add_particle(Particle {
        name: String::from("Sol"),
        radius: 0.0,
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
        radius: 0.0,
        position: Vector3 {
            x: -5.394668400177522,
            y: -7.99940060826745e-1,
            z: 1.24034269478065e-1,
        },
        velocity: Vector3 {
            x: 1.017430886009646e-3 * velocity_scale,
            y: -7.114767972297869e-3 * velocity_scale,
            z: 6.788739409026838e-6 * velocity_scale,
        },
        mass: 9.547919384243222e-4,
    });

    particle_system.add_particle(Particle {
        name: String::from("Saturn"),
        radius: 0.0,
        position: Vector3 {
            x: -2.023594923525499,
            y: -9.836589990338709,
            z: 2.515230776758679e-1,
        },
        velocity: Vector3 {
            x: 5.160849399720572e-3 * velocity_scale,
            y: -1.147589054910672e-3 * velocity_scale,
            z: -1.852057559666422e-4 * velocity_scale,
        },
        mass: 2.858859806661029e-4,
    });

    particle_system.add_particle(Particle {
        name: String::from("Uranus"),
        radius: 0.0,
        position: Vector3 {
            x: 1.838604285828665e+1,
            y: 7.723267785152875,
            z: -2.09383548816917e-1,
        },
        velocity: Vector3 {
            x: -1.54992807736344e-3 * velocity_scale,
            y: 3.435924384648363e-3 * velocity_scale,
            z: 3.287752876026161e-5 * velocity_scale,
        },
        mass: 4.3662440433515637e-5,
    });

    particle_system.add_particle(Particle {
        name: String::from("Neptune"),
        radius: 0.0,
        position: Vector3 {
            x: 2.830602595046843e+1,
            y: -9.782837959471639,
            z: -4.508807468243554e-1,
        },
        velocity: Vector3 {
            x: 1.00677931471006e-3 * velocity_scale,
            y: 2.979403494209059e-3 * velocity_scale,
            z: -8.458453553880726e-5 * velocity_scale,
        },
        mass: 5.151389020535497e-5,
    });

    let dt = 0.593; // 5% of jup period
    let mut simulation = SimulationBuilder::new()
        .with_particle_system(particle_system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .set_time_step(dt)
        .build()
        .expect("simulation built");

    simulation.run_until(simulation_time);
}

mod bench {
    extern crate test;
    use crate::solar_system;
    use test::Bencher;

    /// Length of simulation (currently in years, set by G)
    const SIMULATION_TIME: f64 = 4e4;

    #[bench]
    fn bench_solar_system(b: &mut Bencher) {
        b.iter(|| solar_system(SIMULATION_TIME))
    }
}
