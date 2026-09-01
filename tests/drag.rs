use std::f64::consts::{E, PI};

use nsim::{
    force::{GRAVITY, NewtonianGravity, ScalarDrag},
    integration::Leapfrog,
    math_util::Vector3,
    particle::{Particle, ParticleSystem},
    simulation::Simulation,
};

fn simulate_drag(damping_rate: f64, initial_speed: f64, dt: f64) -> f64 {
    let tau = damping_rate.recip();

    let mut simulation = Simulation::new()
        .add_particle(Particle {
            name: String::from("test"),
            radius: 1.0,
            position: Vector3::default(),
            velocity: Vector3 {
                x: initial_speed,
                ..Default::default()
            },
            mass: 0.0,
        })
        .use_integrator(Leapfrog)
        .add_force(ScalarDrag { damping_rate })
        .set_time_step(dt)
        .build();

    simulation.run_until(tau);

    simulation.particles().state().velocities().x[0]
}

#[test]
fn linear_drag_causes_exponential_velocity_decay() {
    const TOLERANCE: f64 = 1e-4;
    const EXPECTED_ORDER: f64 = 1.0;
    const ORDER_TOLERANCE: f64 = 0.1;

    let initial_speed = 10.0;
    let damping_rate = 0.5;
    let tau = f64::recip(damping_rate);

    let expected_velocity = initial_speed / E;

    let coarse_dt = tau * 1e-2;
    let medium_dt = tau * 1e-3;
    let fine_dt = tau * 1e-4;

    let coarse_velocity = simulate_drag(damping_rate, initial_speed, coarse_dt);

    let medium_velocity = simulate_drag(damping_rate, initial_speed, medium_dt);

    let fine_velocity = simulate_drag(damping_rate, initial_speed, fine_dt);

    let coarse_error = (coarse_velocity - expected_velocity).abs() / expected_velocity;

    let medium_error = (medium_velocity - expected_velocity).abs() / expected_velocity;

    let fine_error = (fine_velocity - expected_velocity).abs() / expected_velocity;

    println!("Coarse dt: {coarse_dt}");
    println!("Coarse velocity: {coarse_velocity}");
    println!("Coarse error: {coarse_error:e}");
    println!();

    println!("Medium dt: {medium_dt}");
    println!("Medium velocity: {medium_velocity}");
    println!("Medium error: {medium_error:e}");
    println!();

    println!("Fine dt: {fine_dt}");
    println!("Fine velocity: {fine_velocity}");
    println!("Fine error: {fine_error:e}");

    // The numerical solution should approach the analytical solution
    // as the timestep is reduced.
    assert!(
        medium_error < coarse_error,
        "Reducing dt from {coarse_dt} to {medium_dt} did not reduce error: \
         {coarse_error:e} -> {medium_error:e}"
    );

    assert!(
        fine_error < medium_error,
        "Reducing dt from {medium_dt} to {fine_dt} did not reduce error: \
         {medium_error:e} -> {fine_error:e}"
    );

    // The finest solution should meet our accuracy requirement.
    assert!(
        fine_error <= TOLERANCE,
        "Fine timestep error {fine_error:e} exceeds tolerance {TOLERANCE:e}"
    );

    let coarse_to_medium_order = (coarse_error / medium_error).log10();

    let medium_to_fine_order = (medium_error / fine_error).log10();

    println!("Coarse -> medium order: {coarse_to_medium_order}");
    println!("Medium -> fine order: {medium_to_fine_order}");

    assert!(
        (coarse_to_medium_order - EXPECTED_ORDER).abs() <= ORDER_TOLERANCE,
        "Expected approximately first-order convergence, got {coarse_to_medium_order}"
    );

    assert!(
        (medium_to_fine_order - EXPECTED_ORDER).abs() <= ORDER_TOLERANCE,
        "Expected approximately first-order convergence, got {medium_to_fine_order}"
    );
}

#[test]
fn gravity_and_drag_inspiral_binary() {
    let mut system = ParticleSystem::new();

    let particle1_mass = 1.0;
    let particle2_mass = 1.0;
    let initial_distance = 1.0;

    let damping_rate = 0.01;

    let orbital_period = 2.0
        * PI
        * f64::sqrt(
            (initial_distance * initial_distance * initial_distance)
                / (GRAVITY * (particle1_mass + particle2_mass)),
        );
    let relative_vel =
        (GRAVITY * (particle1_mass + particle2_mass) * f64::recip(initial_distance)).sqrt();
    let dt = orbital_period * 1e-4;

    system.add_particle(Particle {
        name: String::from("test"),
        radius: 1.0,
        position: Vector3 {
            x: initial_distance * 0.5,
            ..Default::default()
        },
        velocity: Vector3 {
            y: -relative_vel * 0.5,
            ..Default::default()
        },
        mass: 1.0,
    });
    system.add_particle(Particle {
        name: String::from("test"),
        radius: 1.0,
        position: Vector3 {
            x: -initial_distance * 0.5,
            ..Default::default()
        },
        velocity: Vector3 {
            y: relative_vel * 0.5,
            ..Default::default()
        },
        mass: 1.0,
    });

    let mut simulation = Simulation::new()
        .with_particle_system(system)
        .use_integrator(Leapfrog)
        .add_force(NewtonianGravity)
        .add_force(ScalarDrag { damping_rate })
        .set_time_step(dt)
        .set_diagnostic_interval(orbital_period)
        .build();

    simulation.run_steps((orbital_period / dt).round() as usize);

    // particles move to smaller orbit, inner orbits have higher orbital velocity so KE actually goes up here.
    // let initial_kinetic_energy = simulation.diagnostics().kinetic_energy()[0];
    // let final_kinetic_energy = simulation.diagnostics().kinetic_energy()[1];

    // println!("Initial kinetic energy: {initial_kinetic_energy}");
    // println!("Final kinetic energy: {final_kinetic_energy}");
    // println!();

    // assert!(
    //     final_kinetic_energy < initial_kinetic_energy,
    //     "Drag should reduce kinetic energy"
    // );

    let initial_angular_momentum = simulation
        .diagnostics()
        .get_sample(0)
        .angular_momentum()
        .norm();
    let final_angular_momentum = simulation
        .diagnostics()
        .get_sample(1)
        .angular_momentum()
        .norm();

    let initial_total_energy = simulation.diagnostics().get_sample(0).total_energy();
    let final_total_energy = simulation.diagnostics().get_sample(1).total_energy();

    let final_distance = (simulation.particles().state().positions().vector_at(0)
        - simulation.particles().state().positions().vector_at(1))
    .norm();

    println!("Initial angular momentum: {initial_angular_momentum}");
    println!("Final angular momentum: {final_angular_momentum}");
    println!();

    assert!(
        final_angular_momentum < initial_angular_momentum,
        "Drag should reduce angular momentum"
    );

    println!("Initial total energy: {initial_total_energy}");
    println!("Final total energy: {final_total_energy}");
    println!();

    assert!(
        final_total_energy < initial_total_energy,
        "Drag should remove mechanical energy from the system"
    );

    println!("Initial separation: {initial_distance}");
    println!("Final separation: {final_distance}");
    println!();

    assert!(
        final_distance < initial_distance,
        "Gravity and drag together should cause orbital contraction"
    );
}
