mod types;

use rand::{distributions::Distribution, thread_rng, Rng};
use std::{cell::RefCell, rc::Rc};
pub use types::*;

const THETA: f64 = 0.5f64;

pub struct Cell {
    pub position: Position,
    pub size: f64,
}
pub struct BHTree {
    cell: Cell,
    mass: f64,
    center_of_mass: Position,
    num_particles: u32,
    particle: Option<Particle>,
    nw: Option<Rc<RefCell<BHTree>>>,
    ne: Option<Rc<RefCell<BHTree>>>,
    sw: Option<Rc<RefCell<BHTree>>>,
    se: Option<Rc<RefCell<BHTree>>>,
}

impl BHTree {
    pub fn init(cell: Cell) -> Self {
        BHTree {
            cell,
            mass: 0f64,
            center_of_mass: Position(0f64, 0f64),
            num_particles: 0,
            particle: None,
            nw: None,
            ne: None,
            sw: None,
            se: None,
        }
    }

    fn get_children(&self) -> Vec<Rc<RefCell<BHTree>>> {
        let mut children: Vec<Rc<RefCell<BHTree>>> = Vec::with_capacity(4);

        if let Some(child) = &self.nw {
            children.push(Rc::clone(child));
        }
        if let Some(child) = &self.ne {
            children.push(Rc::clone(child));
        }
        if let Some(child) = &self.sw {
            children.push(Rc::clone(child));
        }
        if let Some(child) = &self.se {
            children.push(Rc::clone(child));
        }

        children
    }

    fn get_quadrant(&self, particle: &Particle) -> Quadrant {
        if particle.position.0 < self.cell.position.0 {
            if particle.position.1 < self.cell.position.1 {
                Quadrant::SW
            } else {
                Quadrant::NW
            }
        } else if particle.position.1 < self.cell.position.1 {
            Quadrant::SE
        } else {
            Quadrant::NE
        }
    }

    fn create_and_insert(&mut self, particle: Particle, quadrant: Quadrant) {
        let size = self.cell.size / 2f64;

        let compute_new = move |position: Position| {
            let bht = BHTree::init(Cell { position, size });
            Rc::new(RefCell::new(bht))
        };

        match quadrant {
            Quadrant::NW => {
                let t = self.nw.get_or_insert(compute_new(Position(
                    self.cell.position.0 - size,
                    self.cell.position.1 + size,
                )));
                t.borrow_mut().insert_particle(particle);
            }

            Quadrant::NE => {
                let t = self.ne.get_or_insert(compute_new(Position(
                    self.cell.position.0 + size,
                    self.cell.position.1 + size,
                )));
                t.borrow_mut().insert_particle(particle);
            }
            Quadrant::SW => {
                let t = self.sw.get_or_insert(compute_new(Position(
                    self.cell.position.0 - size,
                    self.cell.position.1 - size,
                )));
                t.borrow_mut().insert_particle(particle);
            }
            Quadrant::SE => {
                let t = self.se.get_or_insert(compute_new(Position(
                    self.cell.position.0 + size,
                    self.cell.position.1 - size,
                )));
                t.borrow_mut().insert_particle(particle);
            }
        };
    }

    pub fn insert_particle(&mut self, particle: Particle) {
        match self.num_particles {
            0 => {
                self.particle = Some(particle);
            }
            1 => {
                let q = self.get_quadrant(&particle);
                self.create_and_insert(particle, q);

                // we've just subdivided, so we've got to place the existing particle as well
                let existing = self.particle.unwrap();
                let q = self.get_quadrant(&existing);
                self.create_and_insert(existing, q);
                self.particle = None
            }
            _ => {
                let q = self.get_quadrant(&particle);
                self.create_and_insert(particle, q);
            }
        }
        self.center_of_mass = self.mass * self.center_of_mass + particle.mass * particle.position;
        self.mass += particle.mass;
        self.center_of_mass = self.center_of_mass / self.mass;
        self.num_particles += 1;
    }

    pub fn calculate_force(&self, target_particle: &Particle) -> f64 {
        let mut force = 0f64;
        if self.num_particles == 1 {
            force = self.particle.unwrap().calculate_force(target_particle)
        } else {
            let r = self.center_of_mass.distance_to(target_particle.position);
            if self.cell.size / r < THETA {
                let com = Particle {
                    position: self.center_of_mass,
                    mass: self.mass,
                };
                force = target_particle.calculate_force(&com);
            } else {
                for child in self.get_children() {
                    force += (*child).borrow().calculate_force(target_particle);
                }
            }
        }
        force
    }
}

pub fn get_distribution<D>(n: usize, distribution: D) -> Vec<Particle>
where
    D: Distribution<f64>,
{
    let mut rng = thread_rng();
    let x: Vec<f64> = (&mut rng).sample_iter(&distribution).take(n).collect();
    let y: Vec<f64> = (&mut rng).sample_iter(&distribution).take(n).collect();
    let v: Vec<Particle> = x
        .iter()
        .zip(y.iter())
        .map(|(x, y)| Particle {
            position: (*x, *y).into(),
            mass: 1.0,
        })
        .collect();
    v
}

#[cfg(test)]
mod tests {
    use rand_distr::Normal;

    use crate::{get_distribution, BHTree, Cell, Particle, Position, Quadrant, THETA};

    #[test]
    fn gets_right_quadrant() {
        let bht = BHTree::init(Cell {
            position: Position(0.0, 0.0),
            size: 1.0,
        });

        let nw = Particle {
            position: Position(-0.5, 0.5),
            mass: 1.0,
        };
        let ne = Particle {
            position: Position(0.5, 0.5),
            mass: 1.0,
        };
        let sw = Particle {
            position: Position(-0.5, -0.5),
            mass: 1.0,
        };
        let se = Particle {
            position: Position(0.5, -0.5),
            mass: 1.0,
        };

        assert_eq!(bht.get_quadrant(&nw), Quadrant::NW);
        assert_eq!(bht.get_quadrant(&ne), Quadrant::NE);
        assert_eq!(bht.get_quadrant(&sw), Quadrant::SW);
        assert_eq!(bht.get_quadrant(&se), Quadrant::SE);
    }

    #[test]
    fn inserts_particles() {
        let mut bht = BHTree::init(Cell {
            position: Position(0.0, 0.0),
            size: 1.0,
        });

        let nw = Particle {
            position: Position(-0.75, 0.75),
            mass: 1.0,
        };

        let se = Particle {
            position: Position(0.5, -0.5),
            mass: 1.0,
        };

        bht.insert_particle(se);

        assert_eq!(bht.num_particles, 1);
        assert_eq!(bht.particle.unwrap(), se);

        bht.insert_particle(nw);

        let nwse = Particle {
            position: Position(-0.25, 0.25),
            mass: 3.0,
        };

        bht.insert_particle(nwse);

        assert_eq!(bht.mass, 5.0);
    }

    #[test]
    fn calculates_force() {
        let distribution = Normal::new(0.0, 0.1).unwrap();
        let particles = get_distribution(100, distribution);

        // Construct tree
        let mut bht = BHTree::init(Cell {
            position: Position(0.0, 0.0),
            size: 1.0,
        });

        for particle in &particles {
            bht.insert_particle(*particle)
        }

        // Compute forces
        // bht.compute_mass();
        for particle in &particles {
            // With Barnes-Hut
            let bh_force = bht.calculate_force(particle);

            // With Brute-Force
            let mut bf_force = 0.0;
            for target in &particles {
                bf_force += particle.calculate_force(target);
            }

            if f64::abs(bh_force - bf_force) / bf_force > 0.1 {
                panic!(
                    "Unreasonably large error deviation: got {}, expected {} (theta {})",
                    bh_force, bf_force, THETA
                );
            }
        }
    }
}
