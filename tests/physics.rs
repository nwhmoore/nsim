use nsim::{
    diagnostics::Diagnostics,
    force::{ForceBuffer, gravity::GRAVITY},
    integration::leapfrog_timestep,
    math_util::{Geometry, vector3::Vector3},
    particle::{Particle, ParticleSystem},
};
use std::f64::consts::PI;

#[test]
fn two_body_conservation() {
    let time_start = 0.0;
    let time_one_period = 2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt();
    // 1% of period
    let time_step = time_one_period * 0.01;

    let mut system = ParticleSystem::new_system();

    system.new_particle(Particle {
        name: String::from("Sol"),
        radius: 1.0,
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
        mass: Some(1.0),
    });

    system.new_particle(Particle {
        name: String::from("Jupiter"),
        radius: 1.0,
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
        mass: Some(1.0),
    });

    // ------------------------------------------------------------------------

    let mut time = time_start;
    let mut force_buffer = ForceBuffer::new(system.particle_count());
    let mut diagnostics = Diagnostics::default();

    let initial_evaluation = force_buffer.compute_accelerations(system.state());
    diagnostics.record(time, system.state(), initial_evaluation.potential_energy);

    while time <= time_step * 1_000.0 {
        time += time_step;
        let force_evaluation = leapfrog_timestep(system.state_mut(), &mut force_buffer, time_step);
        diagnostics.record(time, system.state(), force_evaluation.potential_energy);
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
    let relative_change_in_energy = (final_energy - initial_energy).abs() / initial_energy;
    let epsilon = 1_000_000.0 * f64::EPSILON;
    assert!(
        change_in_linear_momentum <= epsilon,
        "\\delta p: {change_in_linear_momentum} > {epsilon}"
    );
    assert!(
        change_in_angular_momentum <= epsilon,
        "\\delta L: {change_in_angular_momentum} > {epsilon}"
    );
    assert!(
        relative_change_in_energy <= epsilon,
        "\\delta E / E: {relative_change_in_energy} > {epsilon}"
    );
}

#[test]
fn leapfrog_convergence() {
    // ---------------------------  INITIAL PARAMETERS ------------------------

    let one_period = 2.0 * PI * (5.0_f64.powf(3.0) / (GRAVITY * 1.0)).sqrt();

    let mut initial_system = ParticleSystem::new_system();

    initial_system.new_particle(Particle {
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
        mass: Some(1.0),
    });

    initial_system.new_particle(Particle {
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
        mass: None,
    });

    // ------------------------------------------------------------------------

    let time_step = one_period * 0.01;
    let all_time_steps = [
        time_step / 100.0,
        time_step / 200.0,
        time_step / 400.0,
        time_step / 800.0,
        time_step / 1_600.0,
        time_step / 3_200.0,
    ];
    let mut errors = Vec::with_capacity(all_time_steps.len());

    for this_time_step in all_time_steps {
        let mut this_system = initial_system.clone();

        let mut force_buffer = ForceBuffer::new(this_system.particle_count());
        let _initial_computation = force_buffer.compute_accelerations(this_system.state());

        let steps = (one_period / this_time_step).round() as usize;
        for _ in 0..steps {
            let _ = leapfrog_timestep(this_system.state_mut(), &mut force_buffer, this_time_step);
        }

        let error_geometry = Geometry::calculate_geometry(
            Vector3 {
                x: 5.0,
                y: 0.0,
                z: 0.0,
            } - this_system.state().positions().value_at(1),
        );
        let error = error_geometry.dist().abs();
        // println!("error: {error:e}");
        errors.push(error);
        // if time_step_idx > 0 {
        //     let last_error = errors[time_step_idx - 1];
        //     let pval = (last_error / error).log2();
        //     println!("log2 ( {last_error:e} / {error:e} ) = {pval:e}");
        // }
    }

    let log_dt = all_time_steps.iter().map(|dt| dt.ln()).collect::<Vec<_>>();
    let log_err = errors.iter().map(|e| e.ln()).collect::<Vec<_>>();
    let n = log_dt.len() as f64;
    let mean_x = log_dt.iter().sum::<f64>() / n;
    let mean_y = log_err.iter().sum::<f64>() / n;
    let slope = log_dt
        .iter()
        .zip(&log_err)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>()
        / log_dt.iter().map(|x| (x - mean_x).powi(2)).sum::<f64>();

    assert!(
        (slope - 2.0).abs() < 0.15,
        "convergence order {slope} not close to 2"
    )
}
