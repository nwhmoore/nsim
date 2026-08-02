//! Bodies and shared body-position interfaces used by the simulation.

use crate::math::{Acceleration, Position, Velocity};

/// A gravitating body that can exert force on other objects.
///
/// `LargeBody` values have mass and interact with other `LargeBody` values.
/// In the current simulation, they are the only sources of gravity and are
/// measured in solar-mass units.
pub struct LargeBody {
    /// Name of the body, also used as the output filename stem.
    pub name: String,
    /// Position in astronomical units (AU).
    pub pos: Position,
    /// Velocity in astronomical units per year (AU/year).
    pub vel: Velocity,
    /// Acceleration in astronomical units per year squared (AU/year²).
    pub acc: Acceleration,
    /// Mass in solar-mass units.
    pub mass: f64,
}

/// A massless test particle affected by `LargeBody` gravity sources.
///
/// `SmallBody` values do not exert gravity on other bodies in the current
/// restricted N-body model.
pub struct SmallBody {
    /// Name of the body, also used as the output filename stem.
    pub name: String,
    /// Position in astronomical units (AU).
    pub pos: Position,
    /// Velocity in astronomical units per year (AU/year).
    pub vel: Velocity,
    /// Acceleration in astronomical units per year squared (AU/year²).
    pub acc: Acceleration,
}

/// Provides access to an object's position.
pub trait Positioned {
    /// Returns the object's current position.
    fn position(&self) -> &Position;
}

impl Positioned for SmallBody {
    fn position(&self) -> &Position {
        &self.pos
    }
}

impl Positioned for LargeBody {
    fn position(&self) -> &Position {
        &self.pos
    }
}
