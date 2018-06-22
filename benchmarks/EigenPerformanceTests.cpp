#include "benchmark/benchmark.h"
#include <Eigen/Core>

//for now these are just dummy benchmarks


template<class matrix>
inline matrix prod(const matrix& A, const matrix& B){
  matrix C;
  for(int i = 0; i<C.RowsAtCompileTime; i++)
    for(int j = 0; j<C.ColsAtCompileTime; j++)
      C(i,j) = A.row(i)*B.col(j);
  return C;
}

//custom mult
static void BM_CustomMult22f(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix2f rand1 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f rand2 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f rand3 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f rand4 = Eigen::Matrix2f::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(prod(prod(prod(rand1,rand2),rand3),rand4)+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_CustomMult22f);
static void BM_CustomMult44d(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix4d rand1 = Eigen::Matrix4d::Random();
  Eigen::Matrix4d rand2 = Eigen::Matrix4d::Random();
  Eigen::Matrix4d rand3 = Eigen::Matrix4d::Random();
  Eigen::Matrix4d rand4 = Eigen::Matrix4d::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(prod(prod(prod(rand1,rand2),rand3),rand4)+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_CustomMult44d);

//variable size
static void BM_Eigen_MultX22f(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::MatrixXf rand1 = Eigen::MatrixXf::Random(2,2);
  Eigen::MatrixXf rand2 = Eigen::MatrixXf::Random(2,2);
  Eigen::MatrixXf rand3 = Eigen::MatrixXf::Random(2,2);
  Eigen::MatrixXf rand4 = Eigen::MatrixXf::Random(2,2);
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_MultX22f);
static void BM_Eigen_MultX44d(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::MatrixXd rand1 = Eigen::MatrixXd::Random(4,4);
  Eigen::MatrixXd rand2 = Eigen::MatrixXd::Random(4,4);
  Eigen::MatrixXd rand3 = Eigen::MatrixXd::Random(4,4);
  Eigen::MatrixXd rand4 = Eigen::MatrixXd::Random(4,4);
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_MultX44d);


static void BM_Eigen_Mult22f(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix2f rand1 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f rand2 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f rand3 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f rand4 = Eigen::Matrix2f::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_Mult22f);
static void BM_Eigen_Mult22d(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix2d rand1 = Eigen::Matrix2d::Random();
  Eigen::Matrix2d rand2 = Eigen::Matrix2d::Random();
  Eigen::Matrix2d rand3 = Eigen::Matrix2d::Random();
  Eigen::Matrix2d rand4 = Eigen::Matrix2d::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_Mult22d);

static void BM_Eigen_Mult33f(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix3f rand1 = Eigen::Matrix3f::Random();
  Eigen::Matrix3f rand2 = Eigen::Matrix3f::Random();
  Eigen::Matrix3f rand3 = Eigen::Matrix3f::Random();
  Eigen::Matrix3f rand4 = Eigen::Matrix3f::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_Mult33f);
static void BM_Eigen_Mult33d(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix3d rand1 = Eigen::Matrix3d::Random();
  Eigen::Matrix3d rand2 = Eigen::Matrix3d::Random();
  Eigen::Matrix3d rand3 = Eigen::Matrix3d::Random();
  Eigen::Matrix3d rand4 = Eigen::Matrix3d::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_Mult33d);

static void BM_Eigen_Mult44f(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix4f rand1 = Eigen::Matrix4f::Random();
  Eigen::Matrix4f rand2 = Eigen::Matrix4f::Random();
  Eigen::Matrix4f rand3 = Eigen::Matrix4f::Random();
  Eigen::Matrix4f rand4 = Eigen::Matrix4f::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_Mult44f);
static void BM_Eigen_Mult44d(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix4d rand1 = Eigen::Matrix4d::Random();
  Eigen::Matrix4d rand2 = Eigen::Matrix4d::Random();
  Eigen::Matrix4d rand3 = Eigen::Matrix4d::Random();
  Eigen::Matrix4d rand4 = Eigen::Matrix4d::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1*rand2*rand3*rand4+rand4*rand3*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_Mult44d);

static void BM_Eigen_MultMV9(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix<double,9,1> V9d = Eigen::Matrix<double, 9,1>::Random();
  Eigen::Matrix<double, 9,9> M9d = Eigen::Matrix<double, 9,9>::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(V9d.transpose()*M9d*V9d);
    }
  }
BENCHMARK(BM_Eigen_MultMV9);
static void BM_Eigen_MultMM8(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix<double,8,8> rand1 = Eigen::Matrix<double, 8,8>::Random();
  Eigen::Matrix<double, 8,8> rand2 = Eigen::Matrix<double, 8,8>::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(rand1.transpose()*rand2*rand1);
    }
  }
BENCHMARK(BM_Eigen_MultMM8);
static void BM_Eigen_MultBlock(benchmark::State &state) {
  srand((unsigned int) time(0));
  Eigen::Matrix3d rand = Eigen::Matrix3d::Random();
  Eigen::Matrix<double, 9,9> T4d = Eigen::Matrix<double, 9,9>::Random();
  for (auto _ : state) {
    benchmark::DoNotOptimize(T4d.topLeftCorner<3,3>()*T4d.bottomRightCorner<3,3>()*rand);
    }
  }
BENCHMARK(BM_Eigen_MultBlock);


BENCHMARK_MAIN();
