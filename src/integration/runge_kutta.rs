use crate::{
    force::ForceSystem, integration::Integrator, math_util::Vector3Series, particle::ParticleState,
};

/// 4th order integrator, four force evalutaions per time step
#[derive(Clone, Default)]
pub struct RungeKutta4 {
    intermediate: ParticleState,
    total_dpos: Vector3Series,
    total_dvel: Vector3Series,
}

impl Integrator for RungeKutta4 {
    fn initialize(&mut self, state: &ParticleState) {
        let n = state.particle_count();

        self.intermediate.clone_from(state);

        self.total_dpos = Vector3Series::new_with_zeros(n);
        self.total_dvel = Vector3Series::new_with_zeros(n);
    }

    fn evaluate_timestep(&mut self, state: &mut ParticleState, forces: &mut ForceSystem, dt: f64) {
        // k1
        self.apply_first_stage(state, forces, dt / 2.0, dt / 6.0);
        //k2
        self.apply_intermediate_stage(state, forces, dt / 2.0, dt / 3.0);
        //k3
        self.apply_intermediate_stage(state, forces, dt, dt / 3.0);
        //k4 + total accumulation
        self.apply_final_stage(state, forces,  dt / 6.0);

    }

    fn warn() {}
}

impl RungeKutta4 {
    fn apply_first_stage(
        &mut self,
        state: &ParticleState,
        forces: &mut ForceSystem,
        stage_scale: f64,
        total_scale: f64,
    ) {
        forces.evaluate(state);
        let dvel = forces.buffer().accelerations();

        let n = self.intermediate.particle_count();

        let (state_positions, state_velocities) = state.positions_and_velocities();
        let (intermediate_positions, intermediate_velocities) =
            self.intermediate.positions_and_velocities_mut();

        for i in 0..n {
            let vx = state_velocities.x[i];
            let vy = state_velocities.y[i];
            let vz = state_velocities.z[i];

            self.total_dpos.x[i] = vx * total_scale;
            self.total_dpos.y[i] = vy * total_scale;
            self.total_dpos.z[i] = vz * total_scale;

            self.total_dvel.x[i] = dvel.x[i] * total_scale;
            self.total_dvel.y[i] = dvel.y[i] * total_scale;
            self.total_dvel.z[i] = dvel.z[i] * total_scale;

            // Construct the NEXT intermediate state directly from
            // the original state.
            intermediate_positions.x[i] = state_positions.x[i] + vx * stage_scale;
            intermediate_positions.y[i] = state_positions.y[i] + vy * stage_scale;
            intermediate_positions.z[i] = state_positions.z[i] + vz * stage_scale;

            intermediate_velocities.x[i] = state_velocities.x[i] + dvel.x[i] * stage_scale;
            intermediate_velocities.y[i] = state_velocities.y[i] + dvel.y[i] * stage_scale;
            intermediate_velocities.z[i] = state_velocities.z[i] + dvel.z[i] * stage_scale;
        }
    }

    fn apply_intermediate_stage(
        &mut self,
        state: &ParticleState,
        forces: &mut ForceSystem,
        stage_scale: f64,
        total_scale: f64,
    ) {
        forces.evaluate(&self.intermediate);
        let dvel = forces.buffer().accelerations();

        let n = self.intermediate.particle_count();

        let (state_positions, state_velocities) = state.positions_and_velocities();
        let (intermediate_positions, intermediate_velocities) =
            self.intermediate.positions_and_velocities_mut();

        for i in 0..n {
            let vx = intermediate_velocities.x[i];
            let vy = intermediate_velocities.y[i];
            let vz = intermediate_velocities.z[i];

            self.total_dpos.x[i] += vx * total_scale;
            self.total_dpos.y[i] += vy * total_scale;
            self.total_dpos.z[i] += vz * total_scale;

            self.total_dvel.x[i] += dvel.x[i] * total_scale;
            self.total_dvel.y[i] += dvel.y[i] * total_scale;
            self.total_dvel.z[i] += dvel.z[i] * total_scale;

            // Construct the NEXT intermediate state directly from
            // the original state.
            intermediate_positions.x[i] = state_positions.x[i] + vx * stage_scale;
            intermediate_positions.y[i] = state_positions.y[i] + vy * stage_scale;
            intermediate_positions.z[i] = state_positions.z[i] + vz * stage_scale;

            intermediate_velocities.x[i] = state_velocities.x[i] + dvel.x[i] * stage_scale;
            intermediate_velocities.y[i] = state_velocities.y[i] + dvel.y[i] * stage_scale;
            intermediate_velocities.z[i] = state_velocities.z[i] + dvel.z[i] * stage_scale;
        }
    }

    fn apply_final_stage(
    &mut self,
    state: &mut ParticleState,
    forces: &mut ForceSystem,
    total_scale: f64,
) {
    forces.evaluate(&self.intermediate);
    let dvel = forces.buffer().accelerations();
    let n = self.intermediate.particle_count();
    let intermediate_velocities = self.intermediate.velocities();
    let (state_positions, state_velocities) = state.positions_and_velocities_mut();

    for i in 0..n {
        state_positions.x[i] += self.total_dpos.x[i] + intermediate_velocities.x[i] * total_scale;
        state_positions.y[i] += self.total_dpos.y[i] + intermediate_velocities.y[i] * total_scale;
        state_positions.z[i] += self.total_dpos.z[i] + intermediate_velocities.z[i] * total_scale;

        state_velocities.x[i] += self.total_dvel.x[i] + dvel.x[i] * total_scale;
        state_velocities.y[i] += self.total_dvel.y[i] + dvel.y[i] * total_scale;
        state_velocities.z[i] += self.total_dvel.z[i] + dvel.z[i] * total_scale;
    }
}
}
