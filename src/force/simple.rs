use crate::{
    force::{Force, ForceEvaluation},
    math_util::{Vector3, kahan::KahanAccumulator},
    particle::ParticleState,
};

#[derive(Clone)]
pub struct ConstantAccel {
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

/// Applies harmonic oscillator potential a_i = -(k / m_i) (x_i - center)
#[derive(Clone)]
pub struct HarmonicPotential {
    pub k: f64,
    pub center: Vector3,
}

impl Force for HarmonicPotential {
    fn evaluate(&self, particle_state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = particle_state.positions();
        let spring_constant = self.k;

        for (i, &mass) in particle_state.masses().iter().enumerate() {
            let dx = positions.x[i] - self.center.x;
            let dy = positions.y[i] - self.center.y;
            let dz = positions.z[i] - self.center.z;

            debug_assert!(mass > 0.0);
            let scale = -spring_constant / mass;

            output.accelerations.x[i] += dx * scale;
            output.accelerations.y[i] += dy * scale;
            output.accelerations.z[i] += dz * scale;
        }
    }

    fn calculate_potential_energy(&self, state: &ParticleState) -> Option<f64> {
        let positions = state.positions();
        let n = state.particle_count();

        let mut potential_energy = KahanAccumulator::default();

        for i in 0..n {
            let dx = positions.x[i] - self.center.x;
            let dy = positions.y[i] - self.center.y;
            let dz = positions.z[i] - self.center.z;

            let r2 = dx * dx + dy * dy + dz * dz;

            potential_energy.add(0.5 * self.k * r2)
        }

        Some(potential_energy.total())
    }
}
