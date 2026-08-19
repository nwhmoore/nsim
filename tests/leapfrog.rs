use nsim::{
    force::{ForceSystem, GRAVITY, NewtonianGravity},
    integration::leapfrog_timestep,
    math_util::{Geometry, vector3::Vector3},
    particle::{Particle, ParticleSystem},
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

        let mut forces: ForceSystem = ForceSystem::new(this_system.particle_count());
        forces.add_pairwise_force(NewtonianGravity);
        let _initial_computation = forces.evaluate(this_system.state());

        let steps = (one_period / this_time_step).round() as usize;
        for _ in 0..steps {
            let _ = leapfrog_timestep(this_system.state_mut(), &mut forces, this_time_step);
        }

        let error_geometry = Geometry::calculate_geometry(Vector3 {
            x: 5.0 - this_system.state().positions().x[1],
            y: 0.0 - this_system.state().positions().y[1],
            z: 0.0 - this_system.state().positions().z[1],
        });
        let error = error_geometry.dist().abs();
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
