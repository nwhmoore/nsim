use nsim::{
    force::NewtonianGravity,
    integration::Leapfrog,
    math_util::Vector3,
    particle::{Particle, ParticleSystem},
    simulation::Simulation,
};

#[test]
fn massless_particle_receives_gravity_from_massive_particle() {
    let mut system = ParticleSystem::new();

    system.add_particle(Particle {
        name: String::from("test1"),
        radius: 0.0,
        position: Vector3::default(),
        velocity: Vector3::default(),
        mass: 2.0,
    });

    system.add_particle(Particle {
        name: String::from("test2"),
        radius: 0.0,
        position: Vector3 {
            x: 2.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    let simulation = Simulation::new()
        .with_particle_system(system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .build();

    let accelerations = simulation.force_system().buffer().accelerations();
    let test_idx = 1;

    assert!((-0.5 - accelerations.x[test_idx]).abs() < 1e-15);
    assert_eq!(accelerations.y[test_idx], 0.0);
    assert_eq!(accelerations.z[test_idx], 0.0);
}

#[test]
fn massless_particle_does_not_accelerate_massive_particle() {
    let mut system = ParticleSystem::new();

    system.add_particle(Particle {
        name: String::from("massive"),
        radius: 0.0,
        position: Vector3::default(),
        velocity: Vector3::default(),
        mass: 2.0,
    });

    system.add_particle(Particle {
        name: String::from("test"),
        radius: 0.0,
        position: Vector3 {
            x: 2.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    let simulation = Simulation::new()
        .with_particle_system(system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .build();

    let accelerations = simulation.force_system().buffer().accelerations();
    let massive_idx = 0;

    assert_eq!(accelerations.x[massive_idx], 0.0);
    assert_eq!(accelerations.y[massive_idx], 0.0);
    assert_eq!(accelerations.z[massive_idx], 0.0);
}

#[test]
fn massless_particles_do_not_gravitate_each_other() {
    let mut system = ParticleSystem::new();

    system.add_particle(Particle {
        name: String::from("massive"),
        radius: 0.0,
        position: Vector3::default(),
        velocity: Vector3::default(),
        mass: 1.0,
    });

    system.add_particle(Particle {
        name: String::from("test1"),
        radius: 0.0,
        position: Vector3 {
            x: 1.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    system.add_particle(Particle {
        name: String::from("test2"),
        radius: 0.0,
        position: Vector3 {
            y: 1.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    let simulation = Simulation::new()
        .with_particle_system(system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .build();

    let accelerations = simulation.force_system().buffer().accelerations();

    let test1_idx = 1;
    let test2_idx = 2;

    // test1 is accelerated only toward the massive particle at the origin.
    assert!((-1.0 - accelerations.x[test1_idx]).abs() < 1e-15);
    assert_eq!(accelerations.y[test1_idx], 0.0);
    assert_eq!(accelerations.z[test1_idx], 0.0);

    // test2 is accelerated only toward the massive particle at the origin.
    assert_eq!(accelerations.x[test2_idx], 0.0);
    assert!((-1.0 - accelerations.y[test2_idx]).abs() < 1e-15);
    assert_eq!(accelerations.z[test2_idx], 0.0);
}

#[test]
fn massless_particles_do_not_contribute_to_gravitational_potential_energy() {
    let mut system = ParticleSystem::new();

    system.add_particle(Particle {
        name: String::from("massive1"),
        radius: 0.0,
        position: Vector3::default(),
        velocity: Vector3::default(),
        mass: 2.0,
    });

    system.add_particle(Particle {
        name: String::from("massive2"),
        radius: 0.0,
        position: Vector3 {
            x: 2.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 3.0,
    });

    system.add_particle(Particle {
        name: String::from("test"),
        radius: 0.0,
        position: Vector3 {
            x: 10.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    let simulation = Simulation::new()
        .with_particle_system(system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .build();

    let potential_energy = simulation
        .diagnostics()
        .records()
        .first()
        .unwrap()
        .potential_energy();

    // U = -G * m1 * m2 / r
    //   = -1 * 2 * 3 / 2
    //   = -3
    assert!((-3.0 - potential_energy).abs() < 1e-15);
}

#[test]
fn system_of_only_massless_particles_has_no_gravity() {
    let mut system = ParticleSystem::new();

    system.add_particle(Particle {
        name: String::from("test1"),
        radius: 0.0,
        position: Vector3 {
            x: 1.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    system.add_particle(Particle {
        name: String::from("test2"),
        radius: 0.0,
        position: Vector3 {
            y: 1.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    system.add_particle(Particle {
        name: String::from("test3"),
        radius: 0.0,
        position: Vector3 {
            z: 1.0,
            ..Default::default()
        },
        velocity: Vector3::default(),
        mass: 0.0,
    });

    let simulation = Simulation::new()
        .with_particle_system(system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .build();

    let accelerations = simulation.force_system().buffer().accelerations();

    for idx in 0..3 {
        assert_eq!(accelerations.x[idx], 0.0);
        assert_eq!(accelerations.y[idx], 0.0);
        assert_eq!(accelerations.z[idx], 0.0);
    }
}
