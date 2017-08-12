#include <string> 
#include <vector>
#include <array>
#include <fstream>
#include <Eigen/Core>

#include "helpers.h"

#include <iostream>

namespace TC
{

//----------------------------------------------------------------------------
// Node
//----------------------------------------------------------------------------
template <typename _Scalar, int _Dim>
class Node{
public:
  typedef _Scalar Scalar;
  enum {Dim = _Dim};
  typedef Eigen::Matrix<Scalar, Dim, 1> VectorType;
protected:
  VectorType _X;
  VectorType _u;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Node() {_u.setZero();} // <-- maybe I should change it for Node(const vec&u = vec::zero()):_u(u){}
  Node(const VectorType& X):_X(X) { _u.setZero(); } // and here too
  const VectorType& X() const { return _X; }
  const VectorType& u() const { return _u; }
  VectorType x() const {return X()+u();}
  void setX(const VectorType& X) { _X = X; }
  void setu(const VectorType& u) { _u = u; }
  void setx(const VectorType& x) { _u = x - X(); }
  };

//----------------------------------------------------------------------------
// Element
//----------------------------------------------------------------------------
template <typename _Scalar, int _Dim, int _NumNodes>
class Element{
public:
  typedef _Scalar Scalar;
  enum { 
      Dim = _Dim,
      NumNodes = _NumNodes,
      NumDofs = _NumNodes*_Dim
      };
  typedef Node<Scalar, Dim> Node;
  typedef std::array<int, Dim> ConnType; 
  typedef Eigen::Matrix<Scalar, NumDofs, 1> VectorType;
  typedef Eigen::Matrix<Scalar, NumDofs, NumDofs> MatrixType;
protected:
  ConnType _conn;
  const std::vector<Node>& _nodes;
public:
  Element(const std::vector<Node>& nodes): _nodes(nodes){}
  const Node& node(int nodeNumber) const {return _nodes[_conn[nodeNumber]];}
  template <typename T>
  void setConn(const T& conn){
    if(conn.size() < _conn.size()) return;
    for(int i = 0; i<_conn.size(); i++) 
      _conn[i] = conn[i];
    }
  };

//Quadratic 3D element with 20 nodes and reduced integration.
//We can template the material because all benchmarks contain contain 
//  models with only one material. It's not a good idea otherwise.
template <typename MatType>
class C3D20R: Element<typename MatType::Scalar, 3, 20>{
public:
  typedef Element<typename MatType::Scalar, 3, 20> Base;
  typedef typename Base::Scalar Scalar;
  typedef typename Base::Node Node;
  typedef typename Base::VectorType VectorType;
  typedef typename Base::VectorType MatrixType;
protected:
  MatType _mat;
public:
  C3D20R(const MatType& mat, const std::vector<Node>& nodes): 
      Base(nodes), _mat(mat)
    {
    EIGEN_STATIC_ASSERT(MatType::Dim == 3, YOU_MADE_A_PROGRAMMING_MISTAKE);
    }
  VectorType F() const;
  //MatrixType K() const;
  //MatrixType M() const;
  //MatrixType LM() const;
  };

//----------------------------------------------------------------------------
// Finite element model
//----------------------------------------------------------------------------
template <typename _Scalar = double, int _Dim = 3>
class FEM{
public:
  typedef _Scalar Scalar;
  enum {
    Dim = _Dim,
    elDim = _Dim == 3 ? 20 : 8 //cubic elements with interior node
    };
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Dim> NodalMatrix;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, elDim> ConnMatrix; //conn for connectivity
protected:
  NodalMatrix _nodes;
  ConnMatrix _conn;
public:
  FEM(){};
  bool ReadAbaqusInp(const std::string& filename);
  };

//----------------------------------------------------------------------------
// function definitions
//----------------------------------------------------------------------------

template <typename MatType>
typename C3D20R<MatType>::VectorType C3D20R<MatType>::F() const {
  //obtaining uIi. Note that this is the transpose of belytschkos notation
  Eigen::Matrix<Scalar, 8, 3> u;
  for (int i = 0; i< 8; i++)
    u.row(i) = this->node(i).u().transpose();
  //storage for the force
  Eigen::Matrix<Scalar, 8, 3> f = Eigen::Matrix<Scalar, 8, 3>::Zero();


  Eigen::Matrix3d H; //displacement gradient (not the deformation gradient)
  Eigen::Matrix3d E; //strain
  //Eigen::Matrix3d F;
  //Eigen::Matrix3d S;
  //Eigen::Matrix3d P;

  /*for (int i = 0; i<ips.size(); i++)
    {
    auto& dhdX = ip_info.at(i).dhdX; //shorten the name
    H.noalias() = dhdX.transpose()*u; //transpose of the definition in Belytschko
    E.noalias() = .5*(H + H.transpose() + H*H.transpose());
    //Ft.noalias() = H+Eigen::Matrix2d::Identity(); //transponse of the deformation gradient
    //S.noalias() = GetMaterial().stress(E);
    //P.noalias() = S*Ft;

    double k = ip_info.at(i).detj * ips.at(i).weight * GetSection().thickness;
    //Fi.noalias() = k*dhdX*P;

    // Commented code above is the same as the following expression. Eigen library
    //optimizations make run the code faster this way. Also the conditional allows 
    //to avoid having to clear the vector fi before using it.
    Fi.noalias() += k*dhdX*GetMaterial().stress3D(E)*(H + Eigen::Matrix3d::Identity());
    }*/

  }







//This function can read the mesh from an Abaqus input file. It's a very 
//  basic implementation so don't get too crazy about it. It blindly trusts
//  that the connectivity matches the value of FEM::elDim.
//Returns true if succeeded
template <typename _Scalar, int _Dim>
bool FEM<_Scalar, _Dim>::ReadAbaqusInp(const std::string& filename){
  //Setup
  std::vector<std::vector<Scalar>> tempNodes;
  std::vector<std::vector<int>> tempEls;
  std::ifstream inf(filename.c_str());
  if(!inf) return false;
  std::string line;
  std::getline(inf, line);
  //Ignore lines until a node command
  while (!strStartsWith(line, "*Node"))
    std::getline(inf, line);
  //Read the first node and check if file is invalid
  std::getline(inf, line);
  if(strStartsWith(line, "*")) return false;
  //get nodes while line is not a command
  while(!strStartsWith(line, "*"))
    {
    tempNodes.push_back(std::vector<Scalar>()); 
    tokenize(tempNodes.back(), line,1);
    std::getline(inf, line);
    }
  //Make sure it's the right dimension
  for(int i = 0; i<tempNodes.size(); i++)
    if(tempNodes[i].size() != Dim) return false;
  //Ignore lines until an Element command
  while (!strStartsWith(line, "*Element"))
    std::getline(inf, line);
  //Read the first Element and check if file is invalid
  std::getline(inf, line);
  if (strStartsWith(line, "*")) return false;
  //get Elements while line is not a command
  while (!strStartsWith(line, "*"))
    {
    tempEls.push_back(std::vector<int>());
    bool endsWithComma;
    //if the line ends with comma we consider it to continue next line
    do {
      endsWithComma = tokenize(tempEls.back(), line, 1);
      std::getline(inf, line);
      } while (endsWithComma);
    }
  //Checking the dimensions are right and that connectivity values are
  //  consistent with the number of nodes. Note that Abaqus uses a base
  //  index equals to 1.
  for(int i = 0; i<tempEls.size(); i++)
    {
    if (tempEls[i].size() != elDim);
      return false;
    for (int j = 0; j<tempEls[i].size(); j++)
      if(tempEls[i][j] > tempNodes.size())
        return false;
    }
  //Build the matrices. Since this is done after reading all the data
  //  and performing all the desired checks it's more unlikely that 
  //  changes are performed on the model during failure.
  _nodes.resize(tempNodes.size(), Dim);
  _conn.resize(tempEls.size(), elDim);
  for (int i = 0; i<tempNodes.size(); i++)
    for (int j = 0; j<Dim; j++)
      _nodes(i, j) = tempNodes[i][j];
  for (int i = 0; i<tempEls.size(); i++)
    for (int j = 0; j<elDim; j++)
      _conn(i, j) = tempEls[i][j] - 1; //converting to index base 0
  return true;
  }


} //namespace TC
