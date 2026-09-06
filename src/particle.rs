//! Particles in the simulation.

use crate::math_util::{Vector3, Vector3Series};

/// Collection of particles stored as catalog metadata and simulation state.
///
/// The catalog and state use a structure-of-arrays layout. Their vectors must
/// remain aligned: a particle at index `i` has its name and metadata in the
/// catalog's index `i` and its position, velocity, and mass in the state's
/// index `i`. Built systems store massive particles first and massless
/// particles second.
#[derive(Default, Clone)]
pub struct ParticleSystem {
    /// Stable catalog metadata for every particle.
    catalog: ParticleCatalog,
    /// Mutable numerical state used by the integrator.
    state: ParticleState,
}

impl ParticleSystem {
    /// Creates a new [`ParticleSystemBuilder`].
    #[must_use]
    pub fn builder() -> ParticleSystemBuilder {
        ParticleSystemBuilder::default()
    }

    fn with_capacity(n: usize) -> Self {
        ParticleSystem {
            catalog: ParticleCatalog {
                id: Vec::with_capacity(n),
                name: Vec::with_capacity(n),
                radius: Vec::with_capacity(n),
            },
            state: ParticleState {
                masses: Vec::with_capacity(n),
                positions: Vector3Series::with_capacity(n),
                velocities: Vector3Series::with_capacity(n),
                massive_count: 0,
                particle_count: 0,
            },
        }
    }

    /// Returns the stable catalog metadata for the particles in this system.
    #[must_use]
    pub fn catalog(&self) -> &ParticleCatalog {
        &self.catalog
    }

    /// Returns a view of the simulation state for this particle system.
    #[must_use]
    pub fn state(&self) -> &ParticleState {
        &self.state
    }

    /// Returns the mutable simulation state for this particle system.
    pub fn state_mut(&mut self) -> &mut ParticleState {
        &mut self.state
    }

    /// Returns the number of particles currently stored in the system.
    #[must_use]
    pub fn particle_count(&self) -> usize {
        self.state.particle_count
    }
}

/// Persistent metadata associated with each particle.
///
/// Every vector is indexed by the same particle index.
#[derive(Default, Clone)]
pub struct ParticleCatalog {
    /// Stable numeric ID assigned when the particle is added.
    id: Vec<usize>,
    /// Particle names, also used as output filename stems.
    name: Vec<String>,
    /// Particle radii.
    radius: Vec<f64>,
}

/// Time-varying numerical state stored for all particles.
#[derive(Default, Clone)]
pub struct ParticleState {
    /// Particle masses; zero denotes a massless test particle.
    masses: Vec<f64>,
    /// Cartesian positions.
    positions: Vector3Series,
    /// Cartesian velocities.
    velocities: Vector3Series,
    /// Number of massive particles.
    massive_count: usize,
    /// Number of particles.
    particle_count: usize,
}

impl ParticleState {
    /// Returns the number of particles currently represented in the state.
    #[must_use]
    pub fn particle_count(&self) -> usize {
        self.particle_count
    }

    /// Returns the per-particle masses, including zero for massless test
    /// particles.
    #[must_use]
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }

    /// Returns the position series for all particles.
    #[must_use]
    pub fn positions(&self) -> &Vector3Series {
        &self.positions
    }

    /// Returns the mutable position series for all particles.
    pub fn positions_mut(&mut self) -> &mut Vector3Series {
        &mut self.positions
    }

    /// Returns the velocity series for all particles.
    #[must_use]
    pub fn velocities(&self) -> &Vector3Series {
        &self.velocities
    }

    /// Returns the mutable velocity series for all particles.
    pub fn velocities_mut(&mut self) -> &mut Vector3Series {
        &mut self.velocities
    }

    /// Returns the position and velocity series for all particles.
    #[must_use]
    pub fn positions_and_velocities(&self) -> (&Vector3Series, &Vector3Series) {
        (&self.positions, &self.velocities)
    }

    /// Returns mutable position and velocity series for all particles.
    pub fn positions_and_velocities_mut(&mut self) -> (&mut Vector3Series, &mut Vector3Series) {
        (&mut self.positions, &mut self.velocities)
    }

    /// Returns the number of massive particles.
    #[must_use]
    pub fn massive_count(&self) -> usize {
        self.massive_count
    }
}

/// Initial metadata and state used to add one particle to a
/// [`ParticleSystemBuilder`].
#[derive(Clone)]
pub struct Particle {
    /// Name of the particle, also used as the output filename stem.
    pub name: String,
    /// Radius of the particle.
    pub radius: f64,
    /// Initial position `(x, y, z)`.
    pub position: Vector3,
    /// Initial velocity `(u, v, w)`.
    pub velocity: Vector3,
    /// Mass.
    pub mass: f64,
}

/// Builder for a [`ParticleSystem`].
///
/// Particles retain their insertion order within each mass category. When the
/// builder is consumed, all massive particles are placed before all massless
/// particles while retaining their stable catalog IDs.
#[derive(Default, Clone)]
pub struct ParticleSystemBuilder {
    massive_particles: Vec<(usize, Particle)>,
    massless_particles: Vec<(usize, Particle)>,
    next_particle_id: usize,
}

impl ParticleSystemBuilder {
    /// Adds a particle to the system being built.
    pub fn add_particle(&mut self, particle: Particle) {
        if particle.mass == 0.0 {
            self.massless_particles
                .push((self.next_particle_id, particle));
        } else {
            self.massive_particles
                .push((self.next_particle_id, particle));
        }

        self.next_particle_id += 1;
    }

    /// Builds a [`ParticleSystem`] from the builder.
    #[must_use]
    pub fn build(self) -> ParticleSystem {
        let ParticleSystemBuilder {
            massive_particles,
            massless_particles,
            ..
        } = self;
        let massive_count = massive_particles.len();
        let particle_count = massive_count + massless_particles.len();

        let mut system = ParticleSystem::with_capacity(particle_count);
        system.state.massive_count = massive_count;
        system.state.particle_count = particle_count;

        for (id, particle) in massive_particles.into_iter().chain(massless_particles) {
            system.catalog.id.push(id);
            system.catalog.name.push(particle.name);
            system.catalog.radius.push(particle.radius);
            system.state.masses.push(particle.mass);
            system.state.positions.push(particle.position);
            system.state.velocities.push(particle.velocity);
        }

        system
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn particle_with_mass(mass: f64) -> Particle {
        Particle {
            name: String::new(),
            radius: 0.0,
            position: Vector3::default(),
            velocity: Vector3::default(),
            mass,
        }
    }

    #[test]
    fn builder_places_massive_particles_first() {
        let mut builder = ParticleSystem::builder();
        builder.add_particle(particle_with_mass(0.0));
        builder.add_particle(particle_with_mass(2.0));
        builder.add_particle(particle_with_mass(0.0));
        builder.add_particle(particle_with_mass(1.0));

        let system = builder.build();

        assert_eq!(system.particle_count(), 4);
        assert_eq!(system.state().massive_count(), 2);
        assert_eq!(system.state().masses(), &[2.0, 1.0, 0.0, 0.0]);
    }
}
