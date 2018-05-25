#include "benchmark/benchmark.h"
#include "isotropicLinear.h"
#include "anisotropicLinear.h"




static void BM_isotropic_F(benchmark::State &state) {
  TC::IsotropicLinear<double, 3> iso(200e9, .3);
  while (state.KeepRunning()) {
    iso.Stress(Eigen::Matrix3d::Identity());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_isotropic_F);

//~~~~~~~~~~~~~~~~
static void BM_isotropic_K(benchmark::State &state) {
  TC::IsotropicLinear<double, 3> iso(200e9, .3);
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(iso.Stiffness());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_isotropic_K);

//~~~~~~~~~~~~~~~~
static void BM_aniso_F(benchmark::State &state) {
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
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(aniso.Stress(Eigen::Matrix3d::Identity()));
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_aniso_F);

//~~~~~~~~~~~~~~~~
static void BM_aniso_K(benchmark::State &state) {
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
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(aniso.Stiffness());
    }
  }

// Register the function as a benchmark
BENCHMARK(BM_aniso_K);

//~~~~~~~~~~~~~~~~

BENCHMARK_MAIN();
