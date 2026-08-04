//! Particles in the simulation.

/// All particles in the system
pub struct ParticleSystem {
    pub catalog: ParticleCatalog,
    pub state: ParticleState,
    next_particle_id: usize,
}

impl ParticleSystem {
    pub fn new() -> Self {
        ParticleSystem {
            catalog: ParticleCatalog {
                id: Vec::new(),
                name: Vec::new(),
                radius: Vec::new(),
                alive: Vec::new(),
            },
            state: ParticleState {
                mass: Vec::new(),
                x: Vec::new(),
                y: Vec::new(),
                z: Vec::new(),
                vx: Vec::new(),
                vy: Vec::new(),
                vz: Vec::new(),
            },
            next_particle_id: 0,
        }
    }

    pub fn new_particle(&mut self, particle: Particle) {
        self.catalog.id.push(self.next_particle_id);
        self.next_particle_id += 1;

        self.catalog.name.push(particle.name);
        self.catalog.radius.push(particle.radius);
        self.catalog.alive.push(true);

        self.state.mass.push(particle.mass);

        self.state.x.push(particle.pos.0);
        self.state.y.push(particle.pos.1);
        self.state.z.push(particle.pos.2);

        self.state.vx.push(particle.vel.0);
        self.state.vy.push(particle.vel.1);
        self.state.vz.push(particle.vel.2);
    }
}

pub struct ParticleCatalog {
    // ID number
    pub id: Vec<usize>,
    // Particle name
    pub name: Vec<String>,
    // Particle radius.
    pub radius: Vec<f64>,
    // Particle status
    alive: Vec<bool>,
}

pub struct ParticleState {
    // Particle masses. Particles with mass `None` are massless test particles.
    pub mass: Vec<Option<f64>>,

    // Positions
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,

    // Velocities
    pub vx: Vec<f64>,
    pub vy: Vec<f64>,
    pub vz: Vec<f64>,
}

/// Individual particle
pub struct Particle {
    /// Name of the particle, also used as the output filename stem.
    pub name: String,
    /// Radius of particle
    pub radius: f64,
    /// Position
    pub pos: (f64, f64, f64),
    /// Velocity
    pub vel: (f64, f64, f64),
    /// Mass
    pub mass: Option<f64>,
}
