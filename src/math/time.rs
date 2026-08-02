use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]

/// Time in years
pub struct Time(f64);

impl Time {
    pub fn years(value: f64) -> Self {
        Self(value)
    }

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
