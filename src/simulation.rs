use crate::{
    diagnostics::Diagnostics,
    force::{Force, ForceConfiguration, ForceSystem},
    integration::{Integrator, NoIntegrator},
    particle::{Particle, ParticleSystem},
    time::Time,
};

#[derive(Clone)]
pub struct SimulationBuilder<I: Integrator> {
    particles: ParticleSystem,
    time: Time,
    integrator: I,
    force_config: ForceConfiguration,
    diagnostics: Diagnostics,
}

impl<I: Integrator> SimulationBuilder<I> {
    // TODO: Make infallible
    pub fn build(mut self) -> Simulation<I> {
        let particle_count = self.particles.particle_count();

        I::warn();
        self.integrator.initialize(self.particles.state());

        let mut sim = Simulation {
            particles: self.particles,
            time: self.time,
            integrator: self.integrator,
            forces: ForceSystem::new(self.force_config, particle_count),
            diagnostics: self.diagnostics,
        };

        // pre-allocate the force buffer
        sim.forces.evaluate(sim.particles.state());

        // record initial state
        sim.diagnostics.record_current_state(
            sim.time.current_time,
            sim.particles.state(),
            sim.forces.configured_forces(),
        );

        sim
    }

    pub fn with_particle_system(mut self, particle_system: ParticleSystem) -> Self {
        self.particles = particle_system;
        self
    }

    pub fn add_particle(mut self, particle: Particle) -> Self {
        self.particles.add_particle(particle);
        self
    }

    pub fn add_force<F: Force + 'static>(mut self, force: F) -> Self {
        self.force_config.add_force(force);
        self
    }

    pub fn set_diagnostic_interval(mut self, dt: f64) -> Self {
        self.time.set_diagnostic_interval(dt);
        self
    }

    pub fn set_time_step(mut self, dt: f64) -> Self {
        self.time.time_step = dt;
        self
    }
}

impl SimulationBuilder<NoIntegrator> {
    // TODO: examine RK4 api
    pub fn use_integrator<I: Integrator>(self, integrator: I) -> SimulationBuilder<I> {
        SimulationBuilder {
            particles: self.particles,
            time: self.time,
            integrator,
            force_config: self.force_config,
            diagnostics: self.diagnostics,
        }
    }
}

pub struct Simulation<I: Integrator> {
    particles: ParticleSystem,
    time: Time,
    integrator: I,
    forces: ForceSystem,
    diagnostics: Diagnostics,
}

impl Simulation<NoIntegrator> {
    #[allow(clippy::new_ret_no_self)]
    pub fn new() -> SimulationBuilder<NoIntegrator> {
        SimulationBuilder {
            particles: ParticleSystem::default(),
            time: Time::default(),
            integrator: NoIntegrator,
            force_config: ForceConfiguration::default(),
            diagnostics: Diagnostics::default(),
        }
    }
}

impl<I: Integrator> Simulation<I> {
    pub fn run_one_step(&mut self) {
        self.integrator.evaluate_timestep(
            self.particles.state_mut(),
            &mut self.forces,
            self.time.time_step,
        );

        self.time.current_time += self.time.time_step;

        if self.time.current_time >= self.time.diagnostic_schedule.next_diagnostic_record {
            self.diagnostics.record_current_state(
                self.time.current_time,
                self.particles.state(),
                self.forces.configured_forces(),
            );

            self.time.diagnostic_schedule.next_diagnostic_record +=
                self.time.diagnostic_schedule.diagnostic_interval;
        }
    }

    pub fn run_steps(&mut self, steps: usize) {
        for _ in 0..steps {
            self.run_one_step();
        }
    }

    pub fn run_until(&mut self, end_time: f64) {
        while self.current_time() + self.time.time_step <= end_time {
            self.run_one_step();
        }
    }

    pub fn particles(&self) -> &ParticleSystem {
        &self.particles
    }

    pub fn particles_mut(&mut self) -> &mut ParticleSystem {
        &mut self.particles
    }

    pub fn diagnostics(&self) -> &Diagnostics {
        &self.diagnostics
    }

    pub fn get_time_step(&self) -> f64 {
        self.time.time_step
    }

    pub fn current_time(&self) -> f64 {
        self.time.current_time
    }

    pub fn force_system(&self) -> &ForceSystem {
        &self.forces
    }
}

#[cfg(test)]
mod test {
    use crate::{integration::NoIntegrator, simulation::Simulation};

    #[test]
    fn run_until_runs_exact_number_of_complete_steps() {
        let dt = 0.1;

        let mut simulation = Simulation::new()
            .use_integrator(NoIntegrator)
            .set_time_step(dt)
            .build();

        simulation.run_until(1.0);

        assert!(
            (simulation.current_time() - 1.0).abs() < 1e-12,
            "expected time 1.0, got {}",
            simulation.current_time()
        );
    }

    #[test]
    fn run_until_does_not_exceed_end_time() {
        let dt = 0.1;

        let mut simulation = Simulation::new()
            .use_integrator(NoIntegrator)
            .set_time_step(dt)
            .build();

        simulation.run_until(1.05);

        assert!(simulation.current_time() <= 1.05);
        assert!((simulation.current_time() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn run_until_does_not_advance_past_end_time() {
        let dt = 0.1;

        let mut simulation = Simulation::new()
            .use_integrator(NoIntegrator)
            .set_time_step(dt)
            .build();

        simulation.run_until(0.05);

        assert!((simulation.current_time() - 0.0).abs() < 1e-12);
    }
}
