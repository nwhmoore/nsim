use crate::{
    force::{Force, ForceEvaluation},
    math_util::{KahanAccumulator, Vector3},
    particle::ParticleState,
};

/// Applies a harmonic potential with acceleration `-(k / m) * (position -
/// center)`.
#[derive(Clone)]
pub struct HarmonicPotential {
    /// Spring constant.
    pub k: f64,
    /// Center of the potential.
    pub center: Vector3,
}

impl Force for HarmonicPotential {
    fn evaluate(&self, particle_state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = particle_state.positions();
        let spring_constant = self.k;
        // let massive_indices = particle_state.massive_indices();
        let mass = particle_state.masses();
        let massive_count = particle_state.massive_count();

        for i in 0..massive_count {
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
        let massive_count = state.massive_count();

        let mut potential_energy = KahanAccumulator::default();

        for i in 0..massive_count {
            let dx = positions.x[i] - self.center.x;
            let dy = positions.y[i] - self.center.y;
            let dz = positions.z[i] - self.center.z;

            let r2 = dx * dx + dy * dy + dz * dz;

            potential_energy.add(0.5 * self.k * r2);
        }

        Some(potential_energy.total())
    }
}
