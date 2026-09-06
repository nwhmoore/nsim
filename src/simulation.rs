//! This handles simulation construction and running, orchestrating the force
//! system, integrator, particle system, and diagnostics.

use crate::{
    diagnostics::Diagnostics,
    force::{Force, ForceConfiguration, ForceSystem},
    integration::{Integrator, NoIntegrator},
    particle::{Particle, ParticleSystem},
    time::Time,
};

/// Builder for configuring a simulation.
#[derive(Clone)]
pub struct SimulationBuilder<I: Integrator> {
    particles: ParticleSystem,
    time: Time,
    integrator: I,
    force_config: ForceConfiguration,
    diagnostics: Diagnostics,
}

impl<I: Integrator> SimulationBuilder<I> {
    /// Builds the simulation, evaluates its initial forces, and records its
    /// initial diagnostics.
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
            sim.time.current,
            sim.particles.state(),
            sim.forces.configured_forces(),
        );

        sim
    }

    /// Replaces the builder's particle system.
    #[allow(clippy::return_self_not_must_use)]
    pub fn with_particle_system(mut self, particle_system: ParticleSystem) -> Self {
        self.particles = particle_system;
        self
    }

    /// Adds a particle to the builder's particle system.
    #[allow(clippy::return_self_not_must_use)]
    pub fn add_particle(mut self, particle: Particle) -> Self {
        self.particles.add_particle(particle);
        self
    }

    /// Adds a force to the simulation.
    #[allow(clippy::return_self_not_must_use)]
    pub fn add_force<F: Force + 'static>(mut self, force: F) -> Self {
        self.force_config.add_force(force);
        self
    }

    /// Sets the interval between diagnostic records.
    #[allow(clippy::return_self_not_must_use)]
    pub fn set_diagnostic_interval(mut self, dt: f64) -> Self {
        self.time.set_diagnostic_interval(dt);
        self
    }

    /// Sets the simulation timestep.
    #[allow(clippy::return_self_not_must_use)]
    pub fn set_time_step(mut self, dt: f64) -> Self {
        self.time.step = dt;
        self
    }
}

impl SimulationBuilder<NoIntegrator> {
    /// Selects the integrator used by the simulation.
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

/// [`Simulation`] runs and orchestrates the force system, integrator, particle
/// system, and diagnostics.
pub struct Simulation<I: Integrator> {
    particles: ParticleSystem,
    time: Time,
    integrator: I,
    forces: ForceSystem,
    diagnostics: Diagnostics,
}

impl Simulation<NoIntegrator> {
    /// Creates a new [`SimulationBuilder`].
    #[allow(clippy::new_ret_no_self)]
    #[must_use]
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
    /// Advances the simulation by one timestep.
    pub fn run_one_step(&mut self) {
        self.integrator.evaluate_timestep(
            self.particles.state_mut(),
            &mut self.forces,
            self.time.step,
        );

        self.time.current += self.time.step;

        if self.time.current >= self.time.diagnostic_schedule.next_diagnostic_record {
            self.diagnostics.record_current_state(
                self.time.current,
                self.particles.state(),
                self.forces.configured_forces(),
            );

            self.time.diagnostic_schedule.next_diagnostic_record +=
                self.time.diagnostic_schedule.diagnostic_interval;
        }
    }

    /// Advances the simulation by `steps` timesteps.
    pub fn run_steps(&mut self, steps: usize) {
        for _ in 0..steps {
            self.run_one_step();
        }
    }

    /// Advances while another complete timestep would not exceed `end_time`.
    pub fn run_until(&mut self, end_time: f64) {
        while self.current_time() + self.time.step <= end_time {
            self.run_one_step();
        }
    }

    /// Returns the particle system.
    pub fn particles(&self) -> &ParticleSystem {
        &self.particles
    }

    /// Returns mutable access to the particle system.
    pub fn particles_mut(&mut self) -> &mut ParticleSystem {
        &mut self.particles
    }

    /// Returns the recorded diagnostics.
    pub fn diagnostics(&self) -> &Diagnostics {
        &self.diagnostics
    }

    /// Returns the simulation timestep.
    pub fn get_time_step(&self) -> f64 {
        self.time.step
    }

    /// Returns the simulation's current time.
    pub fn current_time(&self) -> f64 {
        self.time.current
    }

    /// Returns the force system.
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
