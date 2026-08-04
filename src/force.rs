use crate::{GRAVITY, particle::ParticleState};

pub fn gravity_acceleration(dimension_dist: f64, dist_squared: f64, attractor_mass: f64) -> f64 {
    // fvec = m1 avec = g m1 m2 / rmag^3 rvec

    -GRAVITY * attractor_mass * dimension_dist / (dist_squared * dist_squared.sqrt())
}

pub struct ForceBuffer {
    ax: Vec<f64>,
    ay: Vec<f64>,
    az: Vec<f64>,
}

impl ForceBuffer{
    pub fn update_accelerations(&mut self, state: ParticleState){
        for (((mass,x),y),z) in state.mass.into_iter().zip(state.x).zip(state.y).zip(state.z){
            
        }
    }
}