use crate::{
    force::{Force, ForceEvaluation},
    math_util::Vector3,
    particle::ParticleState,
};

#[derive(Clone)]
pub struct ConstantAccel {
    pub accel_vec: Vector3,
}

impl Force for ConstantAccel {
    fn evaluate(&self, state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let n = state.particle_count();

        for i in 0..n {
            output.accelerations.x[i] += self.accel_vec.x;
            output.accelerations.y[i] += self.accel_vec.y;
            output.accelerations.z[i] += self.accel_vec.z;
        }
    }
}
