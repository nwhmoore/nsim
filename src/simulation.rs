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

    pub fn new_simulation() -> Self {
        Self {
            particles: ParticleSystem::default(),
            time: Time::default(),
            integrator: None,
            force_config: ForceConfiguration::default(),
            diagnostics: Diagnostics::default(),
        }
    }

    pub fn add_particle(mut self, particle: Particle) -> Self {
        self.particles.new_particle(particle);
        self
    }

    pub fn add_particle_system(mut self, particle_system: ParticleSystem) -> Self {
        // TODO: append instead of truncate if not empty
        self.particles = particle_system;
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

    pub fn set_end_time(mut self, time: f64) -> Self {
        self.time.end_time = time;
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

    pub fn run(&mut self) {
        let steps =
            ((self.time.end_time - self.time.current_time) / self.time.time_step).round() as usize;
        for _ in 0..steps {
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

    pub fn set_end_time(&mut self, time: f64) {
        self.time.end_time = time;
    }

    pub fn current_time(&self) -> f64 {
        self.time.current_time
    }
}
