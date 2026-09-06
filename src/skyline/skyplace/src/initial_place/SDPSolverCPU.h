#ifndef SDP_SOLVER_CPU_H
#define SDP_SOLVER_CPU_H

#include <vector>
#include "cuda_linalg/EigenDef.h"

namespace sdp_solver
{

using namespace cuda_linalg;

class SDPSolverCPU
{
  public:

    SDPSolverCPU(int n);

    void setObjective(EigenSMatrix& C);
    void invertObjectiveSign();
    void addEqualityConstraint(EigenSMatrix& Ai, double bi);
    void addInequalityConstraint(EigenSMatrix& Gi, double hi);

    void solve();

    const EigenVector& getResult() const;

    double getObjectiveValue() const;

  private:

    int n_; // X is n by n matrix
    int m_; // Num Equality   Constraints
    int p_; // Num Inequality Constraints

    double primal_obj_;

    /* ------------ CPU Data ------------ */
    EigenVector solution_;
    std::vector<double> b_;
    std::vector<double> h_;
    EigenSMatrix C_;
    std::vector<EigenSMatrix> A_;
    std::vector<EigenSMatrix> G_;
    /* ---------------------------------- */
};

}

#endif
