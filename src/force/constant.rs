use crate::{
    force::{Force, ForceEvaluation},
    math_util::Vector3,
    particle::ParticleState,
};

/// A force which provides a constant acceleration independent of mass. Used
/// mostly for testing.
#[derive(Clone)]
pub struct ConstantAccel {
    /// the constant acceleration applied
    pub accel_vec: Vector3,
}

impl Force for ConstantAccel {
    fn evaluate(&self, particle_state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let n = particle_state.particle_count();

        for i in 0..n {
            output.accelerations.x[i] += self.accel_vec.x;
            output.accelerations.y[i] += self.accel_vec.y;
            output.accelerations.z[i] += self.accel_vec.z;
        }
    }
}
