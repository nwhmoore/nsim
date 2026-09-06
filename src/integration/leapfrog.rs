use crate::{force::ForceSystem, integration::Integrator, particle::ParticleState};

/// Second-order symplectic integrator using one force evaluation per timestep.
#[derive(Clone)]
pub struct Leapfrog;

impl Integrator for Leapfrog {
    fn initialize(&mut self, _state: &ParticleState) {}

    /// Advances all particles using kick-drift-kick integration. The force
    /// buffer must contain accelerations for the input state and is refreshed
    /// for the output state.
    fn evaluate_timestep(&mut self, state: &mut ParticleState, forces: &mut ForceSystem, dt: f64) {
        let n = state.particle_count();
        let half_dt = 0.5 * dt;

        let accelerations = forces.buffer().accelerations();
        let (positions, velocities) = state.positions_and_velocities_mut();

        for i in 0..n {
            velocities.x[i] += accelerations.x[i] * half_dt;
            velocities.y[i] += accelerations.y[i] * half_dt;
            velocities.z[i] += accelerations.z[i] * half_dt;

            positions.x[i] += velocities.x[i] * dt;
            positions.y[i] += velocities.y[i] * dt;
            positions.z[i] += velocities.z[i] * dt;
        }

        forces.evaluate(state);

        let accelerations = forces.buffer().accelerations();
        let velocities = state.velocities_mut();

        for i in 0..n {
            velocities.x[i] += accelerations.x[i] * half_dt;
            velocities.y[i] += accelerations.y[i] * half_dt;
            velocities.z[i] += accelerations.z[i] * half_dt;
        }
    }

    fn warn() {}
}
