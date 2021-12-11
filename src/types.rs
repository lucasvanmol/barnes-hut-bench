use std::ops::{Div, Mul};

use derive_more::{Add, AddAssign, From, Sub, SubAssign};
use rand::{distributions::Standard, prelude::Distribution};

const G: f64 = 1f64;

#[derive(PartialEq, Debug)]
pub enum Quadrant {
    NW,
    NE,
    SW,
    SE,
}

#[derive(Clone, Copy, Add, AddAssign, Sub, SubAssign, Debug, PartialEq, From)]
pub struct Position(pub f64, pub f64);

impl Distribution<Position> for Standard {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> Position {
        Position(rng.sample(self), rng.sample(self))
    }
}

impl Mul<Position> for f64 {
    type Output = Position;

    fn mul(self, rhs: Position) -> Self::Output {
        Position(rhs.0 * self, rhs.1 * self)
    }
}

impl Div<f64> for Position {
    type Output = Position;

    fn div(self, rhs: f64) -> Self::Output {
        Position(self.0 / rhs, self.1 / rhs)
    }
}

impl Position {
    pub fn distance_to(&self, other: Position) -> f64 {
        f64::sqrt(f64::powi(self.0 - other.0, 2) + f64::powi(self.1 - other.1, 2))
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Particle {
    pub position: Position,
    pub mass: f64,
}

impl Particle {
    pub fn distance_to(&self, other: &Particle) -> f64 {
        self.position.distance_to(other.position)
    }

    pub fn calculate_force(&self, other: &Particle) -> f64 {
        if other != self {
            G * self.mass * other.mass / self.distance_to(other)
        } else {
            0.0
        }
    }
}
