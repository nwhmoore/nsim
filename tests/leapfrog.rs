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

    let ini_time_step = one_period * 0.01;
    let all_time_steps = [
        ini_time_step / 100.0,
        ini_time_step / 200.0,
        ini_time_step / 400.0,
        ini_time_step / 800.0,
        ini_time_step / 1_600.0,
        ini_time_step / 3_200.0,
    ];
    let mut errors = Vec::with_capacity(all_time_steps.len());

    let sim_builder = SimulationBuilder::new_simulation()
        .add_particle_system(initial_system)
        .use_integrator(Leapfrog)
        .add_pairwise_force(NewtonianGravity)
        .set_end_time(one_period)
        .set_diagnostic_interval(one_period);

    for this_time_step in all_time_steps {
        let mut this_simulation = sim_builder
            .clone()
            .set_time_step(this_time_step)
            .build()
            .expect("simulation built");

        this_simulation.run();

        let jup_position = this_simulation.particles().state().positions().value_at(1);

        let error = Vector3 {
            x: 5.0 - jup_position.x,
            y: 0.0 - jup_position.y,
            z: 0.0 - jup_position.z,
        }
        .square()
        .sqrt();
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
