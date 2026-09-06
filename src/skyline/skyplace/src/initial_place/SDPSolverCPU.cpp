#include <cstdio>

#include "SDPSolverCPU.h"
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

std::shared_ptr<ndarray<int,1>>     nint(const std::vector<int>    &X) { return new_array_ptr<int>(X);    }
std::shared_ptr<ndarray<double,1>>  ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }

Matrix::t eigen2fusion(const cuda_linalg::EigenSMatrix& mat)
{
	int N   = mat.rows();
	int NNZ = mat.nonZeros();

	std::vector<int>    rows(NNZ);
	std::vector<int>    cols(NNZ);
	std::vector<double> vals(NNZ);

	int idx = 0;
	for(int i = 0; i < mat.outerSize(); i++)
	{
		for(cuda_linalg::EigenSMatrix::InnerIterator it(mat, i); it; ++it)
		{
			rows[idx] = it.row();
			cols[idx] = it.col();
			vals[idx] = static_cast<double>(it.value());
			idx++;
		}
	}

	return Matrix::sparse(N, N, nint(rows), nint(cols), ndou(vals) );
}

namespace sdp_solver
{

SDPSolverCPU::SDPSolverCPU(int n)
  : n_         (n),
    m_         (0),
    p_         (0),
    primal_obj_(0.0)
{
  C_.resize(n_, n_);
}

void
SDPSolverCPU::setObjective(EigenSMatrix& C)
{
  C_ = C;
}

void
SDPSolverCPU::invertObjectiveSign()
{
  C_ = -1 * C_;
}

void
SDPSolverCPU::addEqualityConstraint(EigenSMatrix& Ai, double bi)
{
  m_++;
  A_.emplace_back(Ai);
  b_.emplace_back(bi);
}

void
SDPSolverCPU::addInequalityConstraint(EigenSMatrix& Gi, double hi)
{
  p_++;
  G_.emplace_back(Gi);
  h_.emplace_back(hi);
}

void
SDPSolverCPU::solve()
{
  assert(A_.size() == m_);
  assert(G_.size() == p_);

	Model::t solver = new Model("SDP Solver"); auto _M = finally([&]() { solver->dispose(); });
  // "finally" <- This is like delete or free of C/C++
  // but I don't know why this has to be called first

  // Setting up the variables
  Variable::t X = solver->variable("X", Domain::inPSDCone(n_));

  // Make Mosek Objective Matrix
  auto objective_mosek = eigen2fusion(C_);

  // Objective
  solver->objective(ObjectiveSense::Minimize, Expr::dot(objective_mosek, X));

  // Equality Constraints
  for(int i = 0; i < m_; i++)
  {
    auto equality_constraint_mosek = eigen2fusion(A_[i]);
    std::string const_name = "Constraint" + std::to_string(i);
    solver->constraint(const_name.c_str(), 
                       Expr::dot(equality_constraint_mosek, X), 
                       Domain::equalsTo(b_[i]));
  }

  solver->setSolverParam("numThreads", "16");
  solver->solve();

  // Get Objective value
  primal_obj_ = solver->primalObjValue();

  // Solution (x vector)
  auto sol = *(X->slice( nint( {0, 1} ), nint( {1, n_} ))->level());

  solution_.resize(n_);
  int idx = 0;
  for(auto& val : sol)
    solution_(idx++) = val;
}

const EigenVector&
SDPSolverCPU::getResult() const
{
  return solution_;
}

double
SDPSolverCPU::getObjectiveValue() const
{
  return primal_obj_;
}

}
