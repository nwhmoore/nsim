use nsim::{
    diagnostics::Diagnostics,
    force::{ForceSystem, GRAVITY, NewtonianGravity},
    integration::leapfrog_timestep,
    math_util::{Geometry, vector3::Vector3},
    particle::{Particle, ParticleSystem},
};
use std::f64::consts::PI;

/// Tests the conservation of linear momentum, angular momentum, and total
/// energy for a two massive particle system in a highly eccentric orbit.
#[test]
fn two_body_conservation() {
    let time_start = 0.0;
    let mut system = ParticleSystem::new_system();

    system.new_particle(Particle {
        name: String::from("Sol"),
        radius: 0.0,
        position: Vector3 {
            x: -0.5,
            y: 0.0,
            z: 0.0,
        },
        velocity: Vector3 {
            x: 0.0,
            y: -0.5,
            z: 0.0,
        },
        mass: 1.0,
    });

    system.new_particle(Particle {
        name: String::from("Jupiter"),
        radius: 0.0,
        position: Vector3 {
            x: 0.5,
            y: 0.0,
            z: 0.0,
        },
        velocity: Vector3 {
            x: 0.0,
            y: 0.5,
            z: 0.0,
        },
        mass: 1.0,
    });

    let pos_geo = Geometry::calculate_geometry(
        system.state().positions().value_at(0) - system.state().positions().value_at(1),
    );
    let vel_vec = system.state().velocities().value_at(0) - system.state().velocities().value_at(1);

    let relative_speed = vel_vec.square().sqrt();
    let relative_position = pos_geo.dist();

    let gravitational_parameter = GRAVITY * 2.0;
    let specific_energy =
        relative_speed * relative_speed * 0.5 - gravitational_parameter / relative_position;

    let semi_major_axis = -gravitational_parameter * 0.5 / specific_energy;
    let one_period = 2.0 * PI * (semi_major_axis.powf(3.0) / gravitational_parameter).sqrt();
    let steps_per_period = 15_000;
    let dt = one_period / steps_per_period as f64;

    // ------------------------------------------------------------------------

    let mut time = time_start;
    let mut forces: ForceSystem = ForceSystem::new(system.particle_count());
    forces.add_pairwise_force(NewtonianGravity);
    let mut diagnostics = Diagnostics::default();

    let initial_evaluation = forces.evaluate(system.state());
    diagnostics.record(
        time,
        system.state(),
        initial_evaluation.potential_energy.total(),
    );

    for _ in 0..steps_per_period {
        time += dt;
        let force_evaluation = leapfrog_timestep(system.state_mut(), &mut forces, dt);
        diagnostics.record(
            time,
            system.state(),
            force_evaluation.potential_energy.total(),
        );
    }

    let initial_linear_momentum = diagnostics.linear_momentum().value_at(0);
    let final_linear_momentum = diagnostics
        .linear_momentum()
        .value_at(diagnostics.number_samples() - 1);
    let change_in_linear_momentum = (final_linear_momentum - initial_linear_momentum)
        .square()
        .sqrt();

    let initial_angular_momentum = diagnostics.angular_momentum().value_at(0);
    let final_angular_momentum = diagnostics
        .angular_momentum()
        .value_at(diagnostics.number_samples() - 1);
    let change_in_angular_momentum = (final_angular_momentum - initial_angular_momentum)
        .square()
        .sqrt();

    let initial_energy = diagnostics.total_energy()[0];
    let final_energy = diagnostics.total_energy()[diagnostics.number_samples() - 1];
    let relative_change_in_energy = (final_energy - initial_energy).abs() / initial_energy.abs();

    // test tolerances
    let linear_momentum_tolerance = 1e-12;
    let angular_momentum_tolerance = 1e-12;
    let energy_tolerance = 1e-6;

    assert!(
        change_in_linear_momentum <= linear_momentum_tolerance,
        "\\delta p: {change_in_linear_momentum:e} > {linear_momentum_tolerance:e}"
    );
    assert!(
        change_in_angular_momentum <= angular_momentum_tolerance,
        "\\delta L: {change_in_angular_momentum:e} > {angular_momentum_tolerance:e}"
    );
    assert!(
        relative_change_in_energy <= energy_tolerance,
        "\\delta E / E: {relative_change_in_energy:e} > {energy_tolerance:e}"
    );
}

/// Following [Chenciner & Montgomery 2000](https://arxiv.org/abs/math/0011268)
/// we recreate the the famous equal-mass three-body periodic orbit.
///
/// The published initial conditions use `G = 1.0`, and have period `T =
/// 6.32591398`. `nsim` uses `G = 4\*PI^2`, so the equivalent velocities are
/// scaled by `2\*PI` and the period is scaled by `1/(2\*PI)`.
#[test]
fn figure_eight_periodic_orbit() {
    // Published period for G = 1.
    let published_period = 6.32591398;

    //nsim gravity
    assert!((GRAVITY - 4.0 * PI * PI).abs() < f64::EPSILON);

    let period = published_period / (2.0 * PI);

    let steps_per_period = 120_000;
    let dt = period / steps_per_period as f64;

    // initial conditions
    let velocity_scale = 2.0 * PI;

    let initial_positions = [
        Vector3 {
            x: 0.97000436,
            y: -0.24308753,
            z: 0.0,
        },
        Vector3 {
            x: -0.97000436,
            y: 0.24308753,
            z: 0.0,
        },
        Vector3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
    ];

    let initial_velocities = [
        Vector3 {
            x: 0.466203685 * velocity_scale,
            y: 0.43236573 * velocity_scale,
            z: 0.0,
        },
        Vector3 {
            x: 0.466203685 * velocity_scale,
            y: 0.43236573 * velocity_scale,
            z: 0.0,
        },
        Vector3 {
            x: -0.93240737 * velocity_scale,
            y: -0.86473146 * velocity_scale,
            z: 0.0,
        },
    ];

    let mut system = ParticleSystem::new_system();
    for i in 0..3 {
        system.new_particle(Particle {
            name: format!("Body{i}"),
            radius: 0.0,
            position: initial_positions[i],
            velocity: initial_velocities[i],
            mass: 1.0,
        });
    }

    let mut forces: ForceSystem = ForceSystem::new(system.particle_count());
    forces.add_pairwise_force(NewtonianGravity);
    let _ = forces.evaluate(system.state());

    // One-third-period choreography test
    let steps_per_third = steps_per_period / 3;

    for _ in 0..steps_per_third {
        leapfrog_timestep(system.state_mut(), &mut forces, dt);
    }

    // test tolerance
    let tolerance = 1.0e-6;

    let permutation = [2usize, 0, 1];
    for (i, &expected_index) in permutation.iter().enumerate() {
        let &position_error = Geometry::calculate_geometry(
            system.state().positions().value_at(i) - initial_positions[expected_index],
        )
        .dist();
        let &velocity_error = Geometry::calculate_geometry(
            system.state().velocities().value_at(i) - initial_velocities[expected_index],
        )
        .dist();

        assert!(
            position_error < tolerance,
            "Body {i} position error after one-third period: \
            {position_error:e} >= {tolerance:e}"
        );

        assert!(
            velocity_error < tolerance,
            "Body {i} velocity error after one-third period: \
            {velocity_error:e} >= {tolerance:e}"
        );
    }

    // Integrate for remaining one period.
    for _ in steps_per_third..steps_per_period {
        leapfrog_timestep(system.state_mut(), &mut forces, dt);
    }

    // assertions
    for i in 0..3 {
        let &position_error = Geometry::calculate_geometry(
            system.state().positions().value_at(i) - initial_positions[i],
        )
        .dist();
        let &velocity_error = Geometry::calculate_geometry(
            system.state().velocities().value_at(i) - initial_velocities[i],
        )
        .dist();

        assert!(
            position_error < tolerance,
            "Body {i} position error after one period: \
            {position_error:e} >= {tolerance:e}"
        );

        assert!(
            velocity_error < tolerance,
            "Body {i} velocity error after one period: \
            {velocity_error:e} >= {tolerance:e}"
        );
    }
}
