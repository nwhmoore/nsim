//! A strongly typed duration measured in years.

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// A duration measured in years.
///
/// The inner floating-point value is private so callers cannot accidentally
/// mix an unlabelled scalar into a time calculation without making the unit
/// conversion explicit.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct Time(f64);

impl Time {
    /// Constructs a duration from a number of years.
    pub fn years(value: f64) -> Self {
        Self(value)
    }

    /// Returns this duration as a floating-point number of years.
    pub fn as_years(self) -> f64 {
        self.0
    }
}

impl Add for Time {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Time {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl AddAssign for Time {
    fn add_assign(&mut self, rhs: Self) {
        self.0 += rhs.0;
    }
}

impl SubAssign for Time {
    fn sub_assign(&mut self, rhs: Self) {
        self.0 -= rhs.0;
    }
}

impl Mul<f64> for Time {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self(self.0 * rhs)
    }
}

impl Div<f64> for Time {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self(self.0 / rhs)
    }
}

impl Neg for Time {
    type Output = Self;

    fn neg(self) -> Self {
        Self(-self.0)
    }
}
