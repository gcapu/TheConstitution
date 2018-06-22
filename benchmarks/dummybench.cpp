#include "benchmark/benchmark.h"
#include <fstream>

#include "isotropicLinear.h"
#include "anisotropicLinear.h"
#include "plasticKinHard2DPE.h"



static void BM_isotropic_F(benchmark::State &state) {
  TC::IsotropicLinear<double, 3> iso(200e9, .3);
  while (state.KeepRunning()) {
    benchmark::DoNotOptimize(iso.Stress(Eigen::Matrix3d::Identity()));
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
static void BM_plasticKinHard_K(benchmark::State &state) {
  std::vector<double> strainList;
  std::vector<double> stressList;
  //using anisotropic material with isotropic matrix
  while (state.KeepRunning()) {
    TC::PlasticKinHard2DPE<double, 3> mat(2, 2, 2, 0.0);
    strainList.clear();
    stressList.clear();
    strainList.push_back(0);
    stressList.push_back(0);
    Eigen::Matrix3d strain;
    for(int i = 0; i<100; i++)
      {
      double dsval = 3.1415/100.*cos(2.*3.1415*i/100.);
      strain<< dsval , 0, 0, 0, 0, 0, 0, 0, 0;
      TC::PlasticKinHard2DPE<double, 3>::ResultType result = mat.Increment(strain);
      strainList.push_back(sin(2.*3.1415*(i+1)/100.)/2.);
      stressList.push_back(result.S()(0,0));
      }
    }
  std::ofstream file("./strain-stress.csv");
  for(int i = 0; i<strainList.size(); i++)
    file << strainList.at(i) << ", " << stressList.at(i) << std::endl;
  }

// Register the function as a benchmark
BENCHMARK(BM_plasticKinHard_K);

//~~~~~~~~~~~~~~~~

BENCHMARK_MAIN();
