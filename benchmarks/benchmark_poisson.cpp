#include "../src/Normal.hpp"
#include "../src/PoissonRecon.hpp"
#include "../src/utils/sampling.hpp"
#include "benchmark/benchmark.h"

constexpr int DEPTH = 8;
constexpr int NUM_POINTS = 1000000;

struct BenchMarkFixture : public benchmark::Fixture {
  NormalApproximations na;
  std::vector<std::array<double, 3>> ball;
  std::vector<std::array<double, 3>> box;
  std::vector<std::array<double, 3>> random;

  void SetUp(const benchmark::State &state) override {
    ball = sample_sphere(NUM_POINTS, 100);
    box = sample_box(100, 150, 75, NUM_POINTS);
    random = rand_points(-100.0, 100.0, NUM_POINTS);
    na = NormalApproximations(random);
  }
};

BENCHMARK_DEFINE_F(BenchMarkFixture, NormalSplatting)(benchmark::State &state) {
  constexpr int DIM = 3;
  int num_points = 10000;

  for (auto _ : state) {
    PoissonRecon poisson(na.vertices(), na.normals(), na.inward_normals());
    poisson.run();
    benchmark::DoNotOptimize(poisson);
  }
}

BENCHMARK_REGISTER_F(BenchMarkFixture, NormalSplatting)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();
