use crate::{
    particle::ParticleState,
    utils::{KahanAccumulator, VectorSeries},
};

#[derive(Default)]
pub struct Diagnostics {
    pub time: Vec<f64>,

    pub total_mass: Vec<f64>,

    pub kinetic_energy: Vec<f64>,
    //pub grav_potential_energy: Vec<f64>,
    //pub total_energy: Vec<f64>,
    pub linear_momentum: VectorSeries,
    pub angular_momentum: VectorSeries,

    pub center_of_mass_position: VectorSeries,
    //pub center_of_mass_velocity: VectorSeries,
}

impl Diagnostics {
    pub fn record(&mut self, time: f64, state: &ParticleState) {
        self.time.push(time);

        let mut total_mass = KahanAccumulator::default();

        let mut kinetic_energy = KahanAccumulator::default();

        let mut momentum_x = KahanAccumulator::default();
        let mut momentum_y = KahanAccumulator::default();
        let mut momentum_z = KahanAccumulator::default();

        let mut angular_momentum_x = KahanAccumulator::default();
        let mut angular_momentum_y = KahanAccumulator::default();
        let mut angular_momentum_z = KahanAccumulator::default();

        let mut mass_position_x = KahanAccumulator::default();
        let mut mass_position_y = KahanAccumulator::default();
        let mut mass_position_z = KahanAccumulator::default();

        for particle_index in 0..state.mass.len() {
            if let Some(mass) = state.mass[particle_index]
                && state.alive[particle_index]
            {
                let x = state.position.x[particle_index];
                let y = state.position.y[particle_index];
                let z = state.position.z[particle_index];

                let vx = state.velocity.x[particle_index];
                let vy = state.velocity.y[particle_index];
                let vz = state.velocity.z[particle_index];

                total_mass.add(mass);
                kinetic_energy.add(0.5 * mass * (vx * vx + vy * vy + vz * vz));

                momentum_x.add(mass * vx);
                momentum_y.add(mass * vy);
                momentum_z.add(mass * vz);

                angular_momentum_x.add(mass * (y * vz - z * vy));
                angular_momentum_y.add(mass * (z * vx - x * vz));
                angular_momentum_z.add(mass * (x * vy - y * vx));

                mass_position_x.add(mass * x);
                mass_position_y.add(mass * y);
                mass_position_z.add(mass * z);
            }
        }

        let mass = total_mass.total();

        self.total_mass.push(mass);
        self.kinetic_energy.push(kinetic_energy.total());

        self.linear_momentum.x.push(momentum_x.total());
        self.linear_momentum.y.push(momentum_y.total());
        self.linear_momentum.z.push(momentum_z.total());

        self.angular_momentum.x.push(angular_momentum_x.total());
        self.angular_momentum.y.push(angular_momentum_y.total());
        self.angular_momentum.z.push(angular_momentum_z.total());

        self.center_of_mass_position
            .x
            .push(mass_position_x.total() / mass);
        self.center_of_mass_position
            .y
            .push(mass_position_y.total() / mass);
        self.center_of_mass_position
            .z
            .push(mass_position_z.total() / mass);
    }
}
