use crate::{
    force::{Force, ForceEvaluation},
    math_util::{KahanAccumulator, Vector3},
    particle::ParticleState,
};

/// Applies harmonic oscillator potential ``a_i`` = -(k / ``m_i``) (``x_i`` - center)
#[derive(Clone)]
pub struct HarmonicPotential {
    /// spring constant
    pub k: f64,
    /// center of potential
    pub center: Vector3,
}

impl Force for HarmonicPotential {
    fn evaluate(&self, particle_state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = particle_state.positions();
        let spring_constant = self.k;
        let massive_indices = particle_state.massive_indices();
        let mass = particle_state.masses();

        for &i in massive_indices {
            let dx = positions.x[i] - self.center.x;
            let dy = positions.y[i] - self.center.y;
            let dz = positions.z[i] - self.center.z;

            debug_assert!(mass[i] > 0.0);
            let scale = -spring_constant / mass[i];

            output.accelerations.x[i] += dx * scale;
            output.accelerations.y[i] += dy * scale;
            output.accelerations.z[i] += dz * scale;
        }
    }

    fn calculate_potential_energy(&self, state: &ParticleState) -> Option<f64> {
        let positions = state.positions();
        let massive_indices = state.massive_indices();

        let mut potential_energy = KahanAccumulator::default();

        for &i in massive_indices {
            let dx = positions.x[i] - self.center.x;
            let dy = positions.y[i] - self.center.y;
            let dz = positions.z[i] - self.center.z;

            let r2 = dx * dx + dy * dy + dz * dz;

            potential_energy.add(0.5 * self.k * r2);
        }

        Some(potential_energy.total())
    }
}
