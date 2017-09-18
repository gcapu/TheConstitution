#include "benchmark/benchmark.h"
#include "minifem.h"
#include "isotropicLinear.h"
#include "anisotropicLinear.h"




static void BM_mini_isotropic_F(benchmark::State &state) {
  mini::FEM<double,3> mymodel;
  TC::IsotropicLinear<double, 3> iso(200e9, .3);
  if(!mymodel.ReadAbaqusInp("../benchmark/input/bar-4x1x1-2el.inp", iso))
    std::cout << "Warning: failed reading input file." << std::endl;
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(mymodel.F());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_mini_isotropic_F);

//~~~~~~~~~~~~~~~~
static void BM_mini_isotropic_K(benchmark::State &state) {
  mini::FEM<double,3> mymodel;
  TC::IsotropicLinear<double, 3> iso(200e9, .3);
  if(!mymodel.ReadAbaqusInp("../benchmark/input/bar-4x1x1-2el.inp", iso))
    std::cout << "Warning: failed reading input file." << std::endl;
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(mymodel.K());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_mini_isotropic_K);

//~~~~~~~~~~~~~~~~
static void BM_mini_aniso_F(benchmark::State &state) {
  mini::FEM<double,3> mymodel;
  Eigen::Matrix<double, 6, 6> K;
  TC::IsotropicLinear<double, 3> li(200e9, .3);
  double lam = li.lambda();
  double mu = li.mu();
  //using anisotropic material with isotropic matrix
  K << lam+2*mu, lam     , lam     , 0,  0,  0,
       lam     , lam+2*mu, lam     , 0,  0,  0,
       lam     , lam     , lam+2*mu, 0,  0,  0,
       0       , 0       , 0       , mu, 0 , 0,
       0       , 0       , 0       , 0 , mu, 0,
       0       , 0       , 0       , 0 , 0 , mu;
  TC::AnisotropicLinear<double, 3> aniso(K);
  if(!mymodel.ReadAbaqusInp("../benchmark/input/bar-4x1x1-2el.inp", aniso))
    std::cout << "Warning: failed reading input file." << std::endl;
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(mymodel.F());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_mini_aniso_F);

//~~~~~~~~~~~~~~~~
static void BM_mini_aniso_K(benchmark::State &state) {
  mini::FEM<double,3> mymodel;
  Eigen::Matrix<double, 6, 6> K;
  TC::IsotropicLinear<double, 3> li(200e9, .3);
  double lam = li.lambda();
  double mu = li.mu();
  //using anisotropic material with isotropic matrix
  K << lam+2*mu, lam     , lam     , 0,  0,  0,
       lam     , lam+2*mu, lam     , 0,  0,  0,
       lam     , lam     , lam+2*mu, 0,  0,  0,
       0       , 0       , 0       , mu, 0 , 0,
       0       , 0       , 0       , 0 , mu, 0,
       0       , 0       , 0       , 0 , 0 , mu;
  TC::AnisotropicLinear<double, 3> aniso(K);
  if(!mymodel.ReadAbaqusInp("../benchmark/input/bar-4x1x1-2el.inp", aniso))
    std::cout << "Warning: failed reading input file." << std::endl;
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(mymodel.K());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_mini_aniso_K);

//~~~~~~~~~~~~~~~~

BENCHMARK_MAIN();
