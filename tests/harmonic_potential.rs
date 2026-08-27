use nsim::{
    force::HarmonicPotential,
    integration::{Leapfrog, NoIntegrator},
    math_util::Vector3,
    particle::Particle,
    simulation::SimulationBuilder,
};
use std::f64::consts::PI;

#[test]
fn harmonic_potential_evaluates_expected_acceleration() {
    let force = HarmonicPotential {
        k: 12.0,
        center: Vector3 {
            x: 1.0,
            y: 2.0,
            z: 3.0,
        },
    };
    let particle = Particle {
        name: String::from("test"),
        radius: 0.0,
        position: Vector3 {
            x: 4.0,
            y: 0.0,
            z: 9.0,
        },
        velocity: Vector3::default(),
        mass: 3.0,
    };

    // test independent of integrator since no steps are run. only initial
    // evaluation is checked.
    let simulation = SimulationBuilder::new()
        .add_particle(particle)
        .add_force(force)
        .use_integrator(NoIntegrator)
        .build()
        .expect("sim built");

    assert_eq!(simulation.force_system().buffer().accelerations.x[0], -12.0);
    assert_eq!(simulation.force_system().buffer().accelerations.y[0], 8.0);
    assert_eq!(simulation.force_system().buffer().accelerations.z[0], -24.0);
}

#[test]
fn harmonic_potential_matches_quarter_period_solution() {
    let position_tolerance = 1e-5;
    let velocity_tolerance = 1e-5;
    let mass = 2.0;
    let spring_constant = 8.0;

    let center = Vector3 {
        x: 10.0,
        y: -3.0,
        z: 7.0,
    };
    let amplitude = 2.0;

    let omega = f64::sqrt(spring_constant / mass);
    let period = 2.0 * PI / omega;

    // Start displaced from the center along x with zero velocity.
    //
    // Exact solution:
    // x(t) = center.x + amplitude * cos(omega * t)
    //
    // At T / 4:
    // position = center
    // velocity_x = -amplitude * omega
    let initial_position = Vector3 {
        x: center.x + amplitude,
        y: center.y,
        z: center.z,
    };

    let initial_velocity = Vector3::default();

    let steps_per_period = 120_000.00;
    let dt = period / steps_per_period;

    let mut simulation = SimulationBuilder::new()
        .add_particle(Particle {
            name: String::from("test"),
            radius: 0.0,
            mass,
            position: initial_position,
            velocity: initial_velocity,
        })
        .add_force(HarmonicPotential {
            k: spring_constant,
            center,
        })
        .set_time_step(dt)
        .use_integrator(Leapfrog)
        .build()
        .expect("simulation build");

    simulation.run_steps((steps_per_period / 4.0) as usize);

    let state = simulation.particles().state();
    let position = state.positions().vector_at(0);
    let velocity = state.velocities().vector_at(0);

    assert!(
        (position.x - center.x).abs() < position_tolerance,
        "x position: expected {}, got {}",
        center.x,
        position.x
    );

    assert!(
        (position.y - center.y).abs() < position_tolerance,
        "y position: expected {}, got {}",
        center.y,
        position.y
    );

    assert!(
        (position.z - center.z).abs() < position_tolerance,
        "z position: expected {}, got {}",
        center.z,
        position.z
    );

    let expected_vx = -amplitude * omega;

    assert!(
        (velocity.x - expected_vx).abs() < velocity_tolerance,
        "x velocity: expected {}, got {}",
        expected_vx,
        velocity.x
    );

    assert!(
        velocity.y.abs() < velocity_tolerance,
        "unexpected y velocity: {}",
        velocity.y
    );

    assert!(
        velocity.z.abs() < velocity_tolerance,
        "unexpected z velocity: {}",
        velocity.z
    );
}

#[test]
fn harmonic_potential_conserves_energy_and_angular_momentum() {
    let mass = 1.0;
    let k = 1.0;
    let center = Vector3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    // For m = k = 1:
    //
    // omega = sqrt(k / m) = 1
    // period = 2π
    //
    // Choose an initial state with both x and y components so that
    // angular momentum is non-zero.
    let initial_position = Vector3 {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let initial_velocity = Vector3 {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };

    let period = 2.0 * PI;
    let dt = period * 0.001;

    let mut simulation = SimulationBuilder::new()
        .add_particle(Particle {
            name: String::from("test"),
            radius: 0.0,
            position: initial_position,
            velocity: initial_velocity,
            mass,
        })
        .add_force(HarmonicPotential { center, k })
        .use_integrator(Leapfrog)
        .set_time_step(dt)
        .build()
        .expect("sim built");

    simulation.run_until(period);

    let initial = simulation.diagnostics().records().first().unwrap();
    let final_record = simulation.diagnostics().records().last().unwrap();

    let initial_energy = initial.total_energy();
    let final_energy = final_record.total_energy();

    let initial_angular_momentum = initial.angular_momentum();
    let final_angular_momentum = final_record.angular_momentum();

    let relative_energy_error = (final_energy - initial_energy).abs() / initial_energy.abs();

    let angular_momentum_error = (final_angular_momentum - initial_angular_momentum).norm();

    assert!(
        relative_energy_error < 1e-6,
        "relative energy error: {relative_energy_error:e}"
    );

    assert!(
        angular_momentum_error < 1e-6,
        "angular momentum error: {angular_momentum_error:e}"
    );
}
