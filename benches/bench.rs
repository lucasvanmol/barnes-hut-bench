use barnes_hut::*;
use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration, Throughput,
};
use rand::distributions::Standard;
use rand_distr::Normal;

fn calculate_forces_bht(particles: &[Particle]) {
    // Construct tree
    let mut bht = BHTree::init(Cell {
        position: Position(0.0, 0.0),
        size: 1.0,
    });

    for particle in particles {
        bht.insert_particle(*particle)
    }

    // Compute forces
    for particle in particles {
        bht.calculate_force(particle);
    }
}

fn calculate_forces_bf(particles: &[Particle]) {
    for particle in particles {
        for target in particles {
            target.calculate_force(particle);
        }
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let nums = [10, 100, 1000, 10000];
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Barnes-Hut");
    group.plot_config(plot_config.clone());

    let distribution = Standard;
    for n in nums {
        let particles = get_distribution(n, distribution);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("bf uniform", n), &n, |b, _| {
            b.iter(|| calculate_forces_bf(black_box(&particles)))
        });
        group.bench_with_input(BenchmarkId::new("bh uniform", n), &n, |b, _| {
            b.iter(|| calculate_forces_bht(black_box(&particles)))
        });
    }

    let distribution = Normal::new(0.0, 0.3).unwrap();
    for n in nums {
        let particles = get_distribution(n, distribution);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("bf normal", n), &n, |b, _| {
            b.iter(|| calculate_forces_bf(black_box(&particles)))
        });
        group.bench_with_input(BenchmarkId::new("bh normal", n), &n, |b, _| {
            b.iter(|| calculate_forces_bht(black_box(&particles)))
        });
    }

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
