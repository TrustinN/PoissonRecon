#include "../src/Emst.hpp"
#include "../src/Normal.hpp"
#include "../src/Octree.hpp"
#include "../src/RiemannianGraph.hpp"
#include "../src/utils/sampling.hpp"
#include "benchmark/benchmark.h"

constexpr int DEPTH = 8;
constexpr int NUM_POINTS = 1000000;

struct BenchMarkFixture : public benchmark::Fixture {
  Octree tree;
  std::vector<std::array<double, 3>> ball;
  std::vector<std::array<double, 3>> box;
  std::vector<std::array<double, 3>> random;

  void SetUp(const benchmark::State &state) override {
    ball = sample_sphere(NUM_POINTS, 100);
    box = sample_box(100, 150, 75, NUM_POINTS);
    random = rand_points(-100.0, 100.0, NUM_POINTS);
    tree = Octree(ball, DEPTH);
  }
};

BENCHMARK_DEFINE_F(BenchMarkFixture, Construct)(benchmark::State &state) {
  constexpr int DIM = 3;
  int num_points = 10000;

  for (auto _ : state) {
    Octree tree(random);
    benchmark::DoNotOptimize(tree);
  }
}

BENCHMARK_DEFINE_F(BenchMarkFixture, kNN)(benchmark::State &state) {
  constexpr int DIM = 3;
  int num_points = 10000;

  for (auto _ : state) {
    for (int i = 0; i < num_points; i++) {
      auto p = ball[i];
      std::vector<int> res = tree.kNearestNeighbors(p, 15);
      benchmark::DoNotOptimize(res);
    }
  }
}

BENCHMARK_DEFINE_F(BenchMarkFixture, kNNGroup)(benchmark::State &state) {
  constexpr int DIM = 3;
  int num_points = 10000;

  for (auto _ : state) {
    std::vector<std::array<double, 3>> queries(ball.begin(),
                                               ball.begin() + num_points);
    std::vector<std::vector<int>> res = tree.kNearestNeighbors(queries, 15);
    benchmark::DoNotOptimize(res);
  }
}

BENCHMARK_DEFINE_F(BenchMarkFixture, kNNGroupFull)(benchmark::State &state) {

  for (auto _ : state) {
    std::vector<std::vector<int>> res = tree.kNearestNeighbors(ball, 15);
    benchmark::DoNotOptimize(res);
  }
}

BENCHMARK_DEFINE_F(BenchMarkFixture, rGConstruct)(benchmark::State &state) {
  constexpr int DIM = 3;
  int num_points = 10000;

  for (auto _ : state) {
    RiemannianGraph rG(ball, tree);
    benchmark::DoNotOptimize(rG);
  }
}

BENCHMARK_DEFINE_F(BenchMarkFixture, EMSTConstruct)(benchmark::State &state) {
  constexpr int DIM = 3;
  int num_points = 10000;

  for (auto _ : state) {
    Emst emst(ball, tree);
    benchmark::DoNotOptimize(emst);
  }
}

BENCHMARK_DEFINE_F(BenchMarkFixture, NormalA)(benchmark::State &state) {
  for (auto _ : state) {
    NormalApproximations na(ball);
    benchmark::DoNotOptimize(na);
  }
}

BENCHMARK_REGISTER_F(BenchMarkFixture, Construct)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(BenchMarkFixture, kNN)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(BenchMarkFixture, kNNGroup)->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(BenchMarkFixture, kNNGroupFull)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(BenchMarkFixture, rGConstruct)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(BenchMarkFixture, EMSTConstruct)
    ->Unit(benchmark::kMillisecond);
BENCHMARK_REGISTER_F(BenchMarkFixture, NormalA)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
