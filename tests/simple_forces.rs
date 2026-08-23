use nsim::{
    force::HarmonicPotential, integration::NoIntegrator, math_util::Vector3, particle::Particle,
    simulation::SimulationBuilder,
};

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
