#include "benchmark/benchmark.h"
#include "minifem.h"
#include "isotropicLinear.h"


static void BM_mini(benchmark::State &state) {
  TC::FEM<double,3> mymodel;
  while (state.KeepRunning()) {
    if(mymodel.ReadAbaqusInp("./benchmark/input/bar-4x1x1-2el.inp"))
      std::cout << "Warning: failed reading input file." << std::endl;
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_mini);

//~~~~~~~~~~~~~~~~
static void BM_isolinear(benchmark::State &state) {
TC::IsotropicLinear<double, 2> mat(200, 0.3);
while (state.KeepRunning()) {
  mat.Stiffness();
  Eigen::Matrix2d E;
  E << 0.1, 0.05, 0.05, 0.05;
  mat.Stress(E);
  }
  }

  // Register the function as a benchmark
  BENCHMARK(BM_isolinear);

  //~~~~~~~~~~~~~~~~
static void BM_isolinear_eval(benchmark::State &state) {
TC::IsotropicLinear<double, 2> mat(200, 0.3);
while (state.KeepRunning()) {
  mat.Stiffness();
  Eigen::Matrix2d E;
  E << 0.1, 0.05, 0.05, 0.05;
  mat.Stress(E).eval();
  }
}

// Register the function as a benchmark
BENCHMARK(BM_isolinear_eval);

//~~~~~~~~~~~~~~~~

// Define another benchmark
static void BM_isolinear_no_opt(benchmark::State &state) {
  TC::IsotropicLinear<double, 2> mat(200, 0.3);
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(mat.Stiffness());
    Eigen::Matrix2d E;
    E << 0.1, 0.05, 0.05, 0.05;
    benchmark::DoNotOptimize(mat.Stress(E).eval());
    }
  }
BENCHMARK(BM_isolinear_no_opt);

BENCHMARK_MAIN();
