use crate::{
    force::ForceSystem, integration::Integrator, math_util::Vector3Series, particle::ParticleState,
};

/// 4th order integrator, four force evalutaions per time step
#[derive(Clone, Default)]
pub struct RungeKutta4 {
    k1: RK4Stage,
    k2: RK4Stage,
    k3: RK4Stage,
    k4: RK4Stage,
    intermediate: ParticleState,
}

impl Integrator for RungeKutta4 {
    fn initialize(&mut self, state: &ParticleState) {
        let n = state.particle_count();

        self.k1 = RK4Stage::new(n);
        self.k2 = RK4Stage::new(n);
        self.k3 = RK4Stage::new(n);
        self.k4 = RK4Stage::new(n);
        self.intermediate = state.clone();
    }

    fn evaluate_timestep(&mut self, state: &mut ParticleState, forces: &mut ForceSystem, dt: f64) {
        // k1 = f(y_n)
        RungeKutta4::evaluate_stage(state, forces, &mut self.k1);

        // k2 = f(y_n + dt/2 * k1)
        self.intermediate.clone_from(state);
        RungeKutta4::apply_stage(&mut self.intermediate, &self.k1, dt * 0.5);

        RungeKutta4::evaluate_stage(&self.intermediate, forces, &mut self.k2);

        // k3 = f(y_n + dt/2 * k2)
        self.intermediate.clone_from(state);
        RungeKutta4::apply_stage(&mut self.intermediate, &self.k2, dt * 0.5);

        RungeKutta4::evaluate_stage(&self.intermediate, forces, &mut self.k3);

        // k4 = f(y_n + dt * k3)
        self.intermediate.clone_from(state);
        RungeKutta4::apply_stage(&mut self.intermediate, &self.k3, dt);

        RungeKutta4::evaluate_stage(&self.intermediate, forces, &mut self.k4);

        // y_{n+1} = y_n + dt/6 * (k1 + 2k2 + 2k3 + k4)
        RungeKutta4::apply_final_update(state, &self.k1, &self.k2, &self.k3, &self.k4, dt);
    }

    fn warn() {}
}

impl RungeKutta4 {
    fn evaluate_stage(state: &ParticleState, forces: &mut ForceSystem, out: &mut RK4Stage) {
        out.dr.clone_from(state.velocities());

        forces.evaluate(state);
        out.dv.clone_from(forces.buffer().accelerations());
    }

    fn apply_stage(intermediate: &mut ParticleState, stage: &RK4Stage, scale: f64) {
        let n = intermediate.particle_count();
        let (positions, velocities) = intermediate.positions_and_velocities_mut();

        for i in 0..n {
            positions.x[i] += stage.dr.x[i] * scale;
            positions.y[i] += stage.dr.y[i] * scale;
            positions.z[i] += stage.dr.z[i] * scale;

            velocities.x[i] += stage.dv.x[i] * scale;
            velocities.y[i] += stage.dv.y[i] * scale;
            velocities.z[i] += stage.dv.z[i] * scale;
        }
    }

    fn apply_final_update(
        state: &mut ParticleState,
        k1: &RK4Stage,
        k2: &RK4Stage,
        k3: &RK4Stage,
        k4: &RK4Stage,
        dt: f64,
    ) {
        let scale = dt / 6.0;
        let n = state.particle_count();
        let (positions, velocities) = state.positions_and_velocities_mut();

        for i in 0..n {
            let dpos_x = (k1.dr.x[i] + k2.dr.x[i] * 2.0 + k3.dr.x[i] * 2.0 + k4.dr.x[i]) * scale;
            let dpos_y = (k1.dr.y[i] + k2.dr.y[i] * 2.0 + k3.dr.y[i] * 2.0 + k4.dr.y[i]) * scale;
            let dpos_z = (k1.dr.z[i] + k2.dr.z[i] * 2.0 + k3.dr.z[i] * 2.0 + k4.dr.z[i]) * scale;

            let dvel_x = (k1.dv.x[i] + k2.dv.x[i] * 2.0 + k3.dv.x[i] * 2.0 + k4.dv.x[i]) * scale;
            let dvel_y = (k1.dv.y[i] + k2.dv.y[i] * 2.0 + k3.dv.y[i] * 2.0 + k4.dv.y[i]) * scale;
            let dvel_z = (k1.dv.z[i] + k2.dv.z[i] * 2.0 + k3.dv.z[i] * 2.0 + k4.dv.z[i]) * scale;

            positions.x[i] += dpos_x;
            positions.y[i] += dpos_y;
            positions.z[i] += dpos_z;

            velocities.x[i] += dvel_x;
            velocities.y[i] += dvel_y;
            velocities.z[i] += dvel_z;
        }
    }
}

#[derive(Clone, Default)]
struct RK4Stage {
    dr: Vector3Series,
    dv: Vector3Series,
}

impl RK4Stage {
    fn new(n: usize) -> Self {
        RK4Stage {
            // TODO: Remove this field to avoid extra allocation/ copying
            dr: Vector3Series::new_with_zeros(n),
            dv: Vector3Series::new_with_zeros(n),
        }
    }
}
