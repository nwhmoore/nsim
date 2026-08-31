//! Coordinate transformations and utilities

use crate::{force::GRAVITY, math_util::Vector3};

/// Keplerian orbital elements
pub struct OrbitalElements {
    /// Semi-major axis
    semi: f64,
    /// Eccentricity
    ecc: f64,
    /// Inclination (radians)
    inc: f64,
    /// Argument of pericenter (radians)
    arg_peri: f64,
    /// Longitude of ascending node (radians)
    long_asc: f64,
    /// Mean anomaly (radians)
    mean_anom: f64,
}

impl OrbitalElements {
    /// Transform keplerian orbital elements to cartesian state vectors
    #[must_use]
    pub fn to_cart(&self, central_mass: f64) -> (Vector3, Vector3) {
        let mu = central_mass * GRAVITY;

        let ecc_anom = solve_keplers_equation(self.ecc, self.mean_anom);
        let (sin_e, cos_e) = ecc_anom.sin_cos();
        let sqrt_one_minus_e2 = (1.0 - self.ecc * self.ecc).sqrt();

        let dist = self.semi * (1.0 - self.ecc * ecc_anom.cos());

        let orbital_pos = Vector3 {
            x: self.semi * (cos_e - self.ecc),
            y: self.semi * sqrt_one_minus_e2 * sin_e,
            z: 0.0,
        };

        let prefac = (mu * self.semi).sqrt() / dist;

        let orbital_vel = Vector3 {
            x: -prefac * sin_e,
            y: prefac * sqrt_one_minus_e2 * cos_e,
            z: 0.0,
        };

        (
            orbital_pos
                .rotate_z(self.arg_peri)
                .rotate_x(self.inc)
                .rotate_z(self.long_asc),
            orbital_vel
                .rotate_z(self.arg_peri)
                .rotate_x(self.inc)
                .rotate_z(self.long_asc),
        )
    }
}

// TODO: make fallible
fn solve_keplers_equation(ecc: f64, mean_anom: f64) -> f64 {
    let mut ecc_anom = mean_anom;
    for _ in 0..100 {
        let f = ecc_anom - ecc * ecc_anom.sin() - mean_anom;
        let fp = 1.0 - ecc * ecc_anom.cos();
        let next_guess = ecc_anom - f / fp;

        if (next_guess - ecc_anom).abs() <= 1e-14 * (1.0 + next_guess.abs()) {
            return next_guess;
        }

        ecc_anom = next_guess;
    }

    ecc_anom
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const TOL: f64 = 1e-12;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < TOL,
            "expected {expected}, got {actual}"
        );
    }

    fn assert_vector_close(actual: Vector3, expected: Vector3) {
        assert!((actual.x - expected.x).abs() < TOL);
        assert!((actual.y - expected.y).abs() < TOL);
        assert!((actual.z - expected.z).abs() < TOL);
    }

    #[test]
    fn circular_orbit_at_zero_anomaly() {
        let elements = OrbitalElements {
            semi: 5.0,
            ecc: 0.0,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: 0.0,
        };

        let (position, velocity) = elements.to_cart(1.0);

        assert_close(position.x, 5.0);
        assert_close(position.y, 0.0);
        assert_close(position.z, 0.0);

        assert_close(velocity.x, 0.0);
        assert_close(velocity.y, (GRAVITY / 5.0).sqrt());
        assert_close(velocity.z, 0.0);
    }

    #[test]
    fn circular_orbit_at_quarter_orbit() {
        let elements = OrbitalElements {
            semi: 5.0,
            ecc: 0.0,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: PI / 2.0,
        };

        let (position, velocity) = elements.to_cart(1.0);

        assert_close(position.x, 0.0);
        assert_close(position.y, 5.0);
        assert_close(position.z, 0.0);

        let speed = (GRAVITY / 5.0).sqrt();

        assert_close(velocity.x, -speed);
        assert_close(velocity.y, 0.0);
        assert_close(velocity.z, 0.0);
    }

    #[test]
    fn circular_orbit_at_half_orbit() {
        let elements = OrbitalElements {
            semi: 5.0,
            ecc: 0.0,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: PI,
        };

        let (position, velocity) = elements.to_cart(1.0);

        assert_close(position.x, -5.0);
        assert_close(position.y, 0.0);
        assert_close(position.z, 0.0);

        assert_close(velocity.x, 0.0);
        assert_close(velocity.y, -((GRAVITY / 5.0).sqrt()));
        assert_close(velocity.z, 0.0);
    }

    #[test]
    fn eccentric_orbit_at_periapsis() {
        let semi = 5.0;
        let ecc = 0.6;
        let central_mass = 1.0;

        let elements = OrbitalElements {
            semi,
            ecc,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: 0.0,
        };

        let (position, velocity) = elements.to_cart(central_mass);

        let expected_r = semi * (1.0 - ecc);
        let mu = central_mass * GRAVITY;
        let expected_v = (mu * (1.0 + ecc) / (semi * (1.0 - ecc))).sqrt();

        assert_close(position.x, expected_r);
        assert_close(position.y, 0.0);
        assert_close(position.z, 0.0);

        assert_close(velocity.x, 0.0);
        assert_close(velocity.y, expected_v);
        assert_close(velocity.z, 0.0);
    }

    #[test]
    fn eccentric_orbit_at_apoapsis() {
        let semi = 5.0;
        let ecc = 0.6;
        let central_mass = 1.0;

        let elements = OrbitalElements {
            semi,
            ecc,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: PI,
        };

        let (position, velocity) = elements.to_cart(central_mass);

        let mu = central_mass * GRAVITY;
        let expected_r = semi * (1.0 + ecc);
        let expected_v = (mu * (1.0 - ecc) / expected_r).sqrt();

        assert_close(position.x, -expected_r);
        assert_close(position.y, 0.0);
        assert_close(position.z, 0.0);

        assert_close(velocity.x, 0.0);
        assert_close(velocity.y, -expected_v);
        assert_close(velocity.z, 0.0);
    }

    #[test]
    fn inclination_rotates_orbit_out_of_xy_plane() {
        let elements = OrbitalElements {
            semi: 5.0,
            ecc: 0.0,
            inc: PI / 2.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: PI / 2.0,
        };

        let (position, velocity) = elements.to_cart(1.0);

        let speed = (GRAVITY / 5.0).sqrt();

        assert_vector_close(
            position,
            Vector3 {
                x: 0.0,
                y: 0.0,
                z: 5.0,
            },
        );

        assert_vector_close(
            velocity,
            Vector3 {
                x: -speed,
                y: 0.0,
                z: 0.0,
            },
        );
    }

    #[test]
    fn orbital_specific_energy_matches_semi_major_axis() {
        let elements = OrbitalElements {
            semi: 5.0,
            ecc: 0.6,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: 1.234,
        };

        let central_mass = 1.0;
        let mu = central_mass * GRAVITY;

        let (position, velocity) = elements.to_cart(central_mass);

        let r2 = position.x * position.x + position.y * position.y + position.z * position.z;

        let v2 = velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z;

        let energy = 0.5 * v2 - mu / r2.sqrt();
        let expected = -mu / (2.0 * elements.semi);

        assert!(
            (energy - expected).abs() < 1e-12,
            "energy = {energy}, expected = {expected}"
        );
    }

    #[test]
    fn angular_momentum_matches_orbital_elements() {
        let elements = OrbitalElements {
            semi: 5.0,
            ecc: 0.6,
            inc: 0.0,
            arg_peri: 0.0,
            long_asc: 0.0,
            mean_anom: 0.8,
        };

        let mu = GRAVITY;

        let (r, v) = elements.to_cart(1.0);

        let h = Vector3 {
            x: r.y * v.z - r.z * v.y,
            y: r.z * v.x - r.x * v.z,
            z: r.x * v.y - r.y * v.x,
        };

        let h_mag = (h.x * h.x + h.y * h.y + h.z * h.z).sqrt();

        let expected = (mu * elements.semi * (1.0 - elements.ecc.powi(2))).sqrt();

        assert!(
            (h_mag - expected).abs() < 1e-12,
            "h = {h_mag}, expected = {expected}"
        );
    }
}
