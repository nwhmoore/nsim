use crate::{
    diagnostics::Diagnostics,
    error::SimError,
    force::{ForceConfiguration, ForceSystem, PairwiseForce},
    integration::{Integrator, Leapfrog},
    particle::{Particle, ParticleSystem},
    time::Time,
};

#[derive(Clone)]
pub struct SimulationBuilder<I: Integrator = Leapfrog> {
    particles: ParticleSystem,
    time: Time,
    integrator: Option<I>,
    force_config: ForceConfiguration,
    diagnostics: Diagnostics,
}

impl Default for SimulationBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl<I: Integrator> SimulationBuilder<I> {
    pub fn build(self) -> Result<Simulation<I>, SimError> {
        let Some(integrator) = self.integrator else {
            return Err(SimError::MissingIntegrator);
        };
        let particle_count = self.particles.particle_count();

        let mut sim = Simulation {
            particles: self.particles,
            time: self.time,
            integrator,
            forces: ForceSystem::new(self.force_config, particle_count),
            diagnostics: self.diagnostics,
        };

        let initial_evaluation = sim.forces.evaluate(sim.particles.state());
        sim.diagnostics.record(
            sim.time.current_time,
            sim.particles.state(),
            initial_evaluation.potential_energy,
        );

        Ok(sim)
    }

    pub fn new() -> Self {
        Self {
            particles: ParticleSystem::default(),
            time: Time::default(),
            integrator: None,
            force_config: ForceConfiguration::default(),
            diagnostics: Diagnostics::default(),
        }
    }

    pub fn with_particle_system(mut self, particle_system: ParticleSystem) -> Self {
        self.particles = particle_system;
        self
    }

    pub fn add_particle(mut self, particle: Particle) -> Self {
        self.particles.new_particle(particle);
        self
    }

    pub fn use_integrator(mut self, integrator: I) -> Self {
        self.integrator = Some(integrator);
        self
    }

    pub fn add_pairwise_force<F: PairwiseForce + 'static>(mut self, force: F) -> Self {
        self.force_config.add_pairwise_force(force);
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

pub struct Simulation<I: Integrator = Leapfrog> {
    particles: ParticleSystem,
    time: Time,
    integrator: I,
    forces: ForceSystem,
    diagnostics: Diagnostics,
}

impl<I: Integrator> Simulation<I> {
    pub fn advance_one_step(&mut self) {
        let force_evaluation = self.integrator.evaluate_timestep(
            self.particles.state_mut(),
            &mut self.forces,
            self.time.time_step,
        );

        self.time.current_time += self.time.time_step;

        if self.time.current_time >= self.time.diagnostic_schedule.next_diagnostic_record {
            self.diagnostics.record(
                self.time.current_time,
                self.particles.state(),
                force_evaluation.potential_energy,
            );

            self.time.diagnostic_schedule.next_diagnostic_record +=
                self.time.diagnostic_schedule.diagnostic_interval;
        }
    }

    pub fn run_steps(&mut self, steps: usize) {
        for _ in 0..steps {
            self.advance_one_step();
        }
    }

    pub fn run_until(&mut self, end_time: f64) {
        while self.current_time() + self.time.time_step <= end_time {
            self.advance_one_step();
        }
    }

    pub fn particles(&self) -> &ParticleSystem {
        &self.particles
    }

    pub fn diagnostics(&self) -> &Diagnostics {
        &self.diagnostics
    }

    pub fn set_time_step(&mut self, dt: f64) {
        self.time.time_step = dt;
    }

    pub fn current_time(&self) -> f64 {
        self.time.current_time
    }
}

#[cfg(test)]
mod test {
    use crate::{integration::NoIntegrator, simulation::SimulationBuilder};

    #[test]
    fn run_until_runs_exact_number_of_complete_steps() {
        let dt = 0.1;

        let mut simulation = SimulationBuilder::new()
            .use_integrator(NoIntegrator)
            .set_time_step(dt)
            .build()
            .expect("simulation built");

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

        let mut simulation = SimulationBuilder::new()
            .use_integrator(NoIntegrator)
            .set_time_step(dt)
            .build()
            .expect("simulation built");

        simulation.run_until(1.05);

        assert!(simulation.current_time() <= 1.05);
        assert!((simulation.current_time() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn run_until_does_not_advance_past_end_time() {
        let dt = 0.1;

        let mut simulation = SimulationBuilder::new()
            .use_integrator(NoIntegrator)
            .set_time_step(dt)
            .build()
            .expect("simulation built");

        simulation.run_until(0.05);

        assert!((simulation.current_time() - 0.0).abs() < 1e-12);
    }
}
