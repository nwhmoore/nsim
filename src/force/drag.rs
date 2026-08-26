use crate::force::Force;

#[derive(Clone)]
pub struct ScalarDrag {
    pub damping_rate: f64,
}

impl Force for ScalarDrag {
    fn evaluate(
        &self,
        state: &crate::particle::ParticleState,
        output: &mut super::ForceEvaluation<'_>,
    ) {
        let velocites = state.velocities();
        let n = state.particle_count();

        for i in 0..n {
            output.accelerations.x[i] += -velocites.x[i] * self.damping_rate;
            output.accelerations.y[i] += -velocites.y[i] * self.damping_rate;
            output.accelerations.z[i] += -velocites.z[i] * self.damping_rate;
        }
    }
}
