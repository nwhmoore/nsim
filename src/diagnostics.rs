//! Whole-system self recorded as structure-of-arrays time series.
//!
//! Global quantities include only active particles with a `Some` mass value.
//! Massless test particles are intentionally excluded from total mass,
//! kinetic energy, momentum, angular momentum, and center-of-mass quantities.

use crate::{
    force::Force,
    math_util::{
        Vector3,
        kahan::{Kahan3, KahanAccumulator},
    },
    particle::ParticleState,
};

/// Time series of global quantities derived from simulation states.
///
/// Every field has one entry per call to [`Diagnostics::record`]. The entry at
/// a given index therefore refers to the time at the same index in [`Self::time`].
#[derive(Default, Clone)]
pub struct Diagnostics {
    records: Vec<DiagnosticRecord>,
}

impl Diagnostics {
    pub fn evaluate_now(
        current_time: f64,
        particle_state: &ParticleState,
        forces: &[Box<dyn Force>],
    ) -> DiagnosticRecord {
        let mut total_mass = KahanAccumulator::default();
        let mut kinetic_energy = KahanAccumulator::default();
        let mut potential_energy = KahanAccumulator::default();
        let mut linear_momentum = Kahan3::default();
        let mut angular_momentum = Kahan3::default();
        let mut mass_position = Kahan3::default();

        let positions = particle_state.positions();
        let velocities = particle_state.velocities();
        let masses = particle_state.masses();
        let n = particle_state.particle_count();

        for i in 0..n {
            let rx = positions.x[i];
            let ry = positions.y[i];
            let rz = positions.z[i];

            let vx = velocities.x[i];
            let vy = velocities.y[i];
            let vz = velocities.z[i];

            let mass = masses[i];
            total_mass.add(mass);

            kinetic_energy.add(0.5 * mass * (vx * vx + vy * vy + vz * vz));

            linear_momentum.add(vx * mass, vy * mass, vz * mass);

            // angular_momentum.add(&(position.cross(&velocity) * mass));
            angular_momentum.add(
                mass * (ry * vz - rz * vy),
                mass * (rz * vx - rx * vz),
                mass * (rx * vy - ry * vx),
            );

            mass_position.add(mass * rx, mass * ry, mass * rz);
        }

        for force in forces {
            if let Some(energy) = force.calculate_potential_energy(particle_state) {
                potential_energy.add(energy);
            }
        }

        let total_mass = total_mass.total();
        let kinetic_energy = kinetic_energy.total();
        let potential_energy = potential_energy.total();
        let total_energy = kinetic_energy + potential_energy;
        let linear_momentum = linear_momentum.total();
        let angular_momentum = angular_momentum.total();
        let center_of_mass_position = mass_position.total() / total_mass;
        let center_of_mass_velocity = linear_momentum / total_mass;

        DiagnosticRecord {
            current_time,
            total_mass,
            kinetic_energy,
            potential_energy,
            total_energy,
            linear_momentum,
            angular_momentum,
            center_of_mass_position,
            center_of_mass_velocity,
        }
    }

    /// Records one diagnostic sample for a simulation state.
    pub fn record_current_state(
        &mut self,
        current_time: f64,
        particle_state: &ParticleState,
        forces: &[Box<dyn Force>],
    ) {
        let diagnostic_record = Diagnostics::evaluate_now(current_time, particle_state, forces);
        self.records.push(diagnostic_record);
    }

    pub fn records(&self) -> &[DiagnosticRecord] {
        &self.records
    }

    pub fn get_sample(&self, idx: usize) -> &DiagnosticRecord {
        &self.records[idx]
    }

    pub fn number_samples(&self) -> usize {
        self.records.len()
    }
}

#[derive(Clone)]
pub struct DiagnosticRecord {
    current_time: f64,
    /// Total mass of active massive bodies.
    total_mass: f64,
    /// Total kinetic energy of active massive bodies.
    kinetic_energy: f64,
    /// Pairwise Newtonian gravitational potential energy.
    potential_energy: f64,
    /// Sum of kinetic and gravitational potential energy.
    total_energy: f64,
    /// Total linear momentum, stored as parallel component series.
    linear_momentum: Vector3,
    /// Total angular momentum about the simulation origin, stored as parallel
    /// component series.
    angular_momentum: Vector3,
    /// Center-of-mass position of the active massive bodies.
    center_of_mass_position: Vector3,
    /// Center-of-mass velocity of the active massive bodies.
    center_of_mass_velocity: Vector3,
}

impl DiagnosticRecord {
    pub fn current_time(&self) -> f64 {
        self.current_time
    }

    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.kinetic_energy
    }

    pub fn potential_energy(&self) -> f64 {
        self.potential_energy
    }

    pub fn total_energy(&self) -> f64 {
        self.total_energy
    }

    pub fn linear_momentum(&self) -> Vector3 {
        self.linear_momentum
    }

    pub fn angular_momentum(&self) -> Vector3 {
        self.angular_momentum
    }

    pub fn center_of_mass_position(&self) -> Vector3 {
        self.center_of_mass_position
    }

    pub fn center_of_mass_velocity(&self) -> Vector3 {
        self.center_of_mass_velocity
    }
}
