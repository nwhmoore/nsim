//! Physical quantities and gravitational calculations used by the simulation.

use crate::{
    body::{LargeBody, Positioned},
    math::time::Time,
};
use std::{
    f64::consts::PI,
    ops::{AddAssign, Mul, Sub},
};

pub mod time;

/// Newton's gravitational constant
///
/// Expressed in the simulation's unit system as AU³ / (solar mass · year²).
/// Its value, `4π²`, makes a one-solar-mass circular orbit with a semi-major
/// axis of one AU have a period of one year.
pub const GRAVITY: f64 = 4.0 * PI * PI;

/// A three-dimensional position measured in astronomical units (AU).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Position {
    /// X coordinate in AU.
    pub x: f64,
    /// Y coordinate in AU.
    pub y: f64,
    /// Z coordinate in AU.
    pub z: f64,
}

impl AddAssign for Position {
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl Sub for Position {
    type Output = Position;

    fn sub(self, rhs: Position) -> Position {
        Position {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Position {
    /// Returns the squared Euclidean magnitude of this position vector.
    ///
    /// This avoids a square root and is useful when only distances need to be
    /// compared or when calculating inverse-cube gravitational factors.
    pub fn magnitude_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
}

/// A three-dimensional velocity measured in AU/year.
#[derive(Debug, Clone, Copy, Default)]
pub struct Velocity {
    /// X component in AU/year.
    pub u: f64,
    /// Y component in AU/year.
    pub v: f64,
    /// Z component in AU/year.
    pub w: f64,
}

impl AddAssign for Velocity {
    fn add_assign(&mut self, rhs: Self) {
        self.u += rhs.u;
        self.v += rhs.v;
        self.w += rhs.w;
    }
}

impl Mul<Time> for Velocity {
    type Output = Position;

    fn mul(self, dt: Time) -> Self::Output {
        Position {
            x: self.u * dt.as_years(),
            y: self.v * dt.as_years(),
            z: self.w * dt.as_years(),
        }
    }
}

/// A three-dimensional acceleration measured in AU/year².
#[derive(Debug, Clone, Copy, Default)]
pub struct Acceleration {
    /// X component in AU/year².
    pub r: f64,
    /// Y component in AU/year².
    pub s: f64,
    /// Z component in AU/year².
    pub t: f64,
}

impl AddAssign for Acceleration {
    fn add_assign(&mut self, rhs: Self) {
        self.r += rhs.r;
        self.s += rhs.s;
        self.t += rhs.t;
    }
}

impl Mul<Time> for Acceleration {
    type Output = Velocity;

    fn mul(self, dt: Time) -> Velocity {
        Velocity {
            u: self.r * dt.as_years(),
            v: self.s * dt.as_years(),
            w: self.t * dt.as_years(),
        }
    }
}

// /// Keplerian orbital elements of an object in relation to a central mass.
// pub struct OrbitalElements {
//     /// Semi-major axis of an orbit in units of AU
//     semi: f64,
//     /// Eccentricity of an orbit
//     ecc: f64,
//     /// Inclination of an orbit relative to the reference plane in units of radians
//     inc: f64,
//     /// Longitude of ascending node in units of radians
//     long_asc_node: f64,
//     /// Argument of pericenter in units of radians
//     arg_peri: f64,
//     /// Mean anomaly of an orbit in units of radians
//     mean_anom: f64,
// }

/// Calculates the acceleration caused by one [`LargeBody`].
///
/// The acceleration points from `body` toward `large` and has magnitude
/// `GRAVITY * large.mass / distance³`. The source and target must not occupy
/// the same position because the inverse-cube factor is singular there.
///
/// This function is generic over any type implementing [`Positioned`], which
/// allows it to operate on either a [`SmallBody`] or a
/// [`LargeBody`].
///
/// [`SmallBody`]: crate::body::SmallBody
pub fn gravity_acceleration<B>(body: &B, large: &LargeBody) -> Acceleration
where
    B: Positioned,
{
    // fvec = m1 avec = g m1 m2 / rmag^3 rvec

    let r_vec = *body.position() - *large.position();
    let r_squared = r_vec.magnitude_squared();
    let accel_prefac = -GRAVITY * large.mass / (r_squared * r_squared.sqrt());

    Acceleration {
        r: r_vec.x * accel_prefac,
        s: r_vec.y * accel_prefac,
        t: r_vec.z * accel_prefac,
    }
}
