use nsim::{
    force::{GRAVITY, NewtonianGravity},
    integration::Leapfrog,
    math_util::Vector3,
    particle::{Particle, ParticleSystem},
    simulation::SimulationBuilder,
};
use std::f64::consts::PI;

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
        mass: 1.0,
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
        mass: 0.0,
    });

    // ------------------------------------------------------------------------

    let steps_per_periods = [10_000, 20_000, 40_000, 80_000, 160_000, 320_000];

    let all_time_steps = steps_per_periods
        .iter()
        .map(|&steps| one_period / steps as f64)
        .collect::<Vec<_>>();

    let mut errors = Vec::with_capacity(steps_per_periods.len());

    let sim_builder = SimulationBuilder::new()
        .with_particle_system(initial_system)
        .use_integrator(Leapfrog)
        .add_pairwise_force(NewtonianGravity)
        .set_diagnostic_interval(one_period);

    for (&step_per_period, &this_time_step) in steps_per_periods.iter().zip(all_time_steps.iter()) {
        let mut this_simulation = sim_builder
            .clone()
            .set_time_step(this_time_step)
            .build()
            .expect("simulation built");

        this_simulation.run_steps(step_per_period);

        let jup_position = this_simulation.particles().state().positions().value_at(1);

        let error = Vector3 {
            x: 5.0 - jup_position.x,
            y: -jup_position.y,
            z: -jup_position.z,
        }
        .norm();
        println!("error: {error:e}");
        errors.push(error);
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
