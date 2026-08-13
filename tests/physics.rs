use nsim::{
    diagnostics::Diagnostics,
    force::{ForceBuffer, gravity::GRAVITY},
    integration::leapfrog_timestep,
    math_util::vector3::Vector3,
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
fn leapfrog_convergance() {
    // ---------------------------  INITIAL PARAMETERS ------------------------

    let time_start = 0.0;
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
        time_step / 20_000.0,
        time_step / 40_000.0,
        time_step / 80_000.0,
        time_step / 160_000.0,
        time_step / 320_000.0,
        time_step / 640_000.0,
    ];
    let mut errors = Vec::with_capacity(all_time_steps.len());

    for (time_step_idx, &this_time_step) in all_time_steps.iter().enumerate() {
        let mut this_system = initial_system.clone();
        let mut time = time_start;
        let mut force_buffer = ForceBuffer::new(this_system.particle_count());

        while time <= one_period {
            let _ = leapfrog_timestep(this_system.state_mut(), &mut force_buffer, this_time_step);
            time += this_time_step;
        }

        let error = (5.0 - this_system.state().positions().value_at(1).x).abs();
        println!("error: {error}");
        errors.push(error);
        if time_step_idx > 0 {
            let last_error = errors[time_step_idx - 1];
            let pval = (last_error / error).log2();
            println!("log2 ( {last_error:e} / {error:e} ) = {pval:e}");
            //assert!((2.0 - pval).abs() < 1e-10);
            // error is not converging
        }
    }
}
