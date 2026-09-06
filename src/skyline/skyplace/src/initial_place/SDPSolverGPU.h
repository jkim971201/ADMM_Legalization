#ifndef SDP_SOLVER_GPU_H
#define SDP_SOLVER_GPU_H

#include <memory>
#include <vector>
#include "cuda_linalg/EigenDef.h"

#include "cuda_linalg/CudaVector.h"
#include "cuda_linalg/CudaMatrix.h"
#include "cuda_linalg/CudaSparseMatrix.h"

namespace sdp_solver
{

using namespace cuda_linalg;

struct LbfgsNode;

struct Parameters
{
  int    max_alm_iter;
  int    max_lbfgs_iter;
  int    log_freq;
  int    lbfgs_length;
  double tol_alm;
  double tol_lbfgs;
  double rho_certificate;
  double rho_update_factor;
  bool   verbose;
};

class SDPSolverGPU
{
  public:

    SDPSolverGPU(int n);

    void setObjective(EigenSMatrix& C);
    void setMaximization();
    void setVerbose(bool is_verbose);
    void addEqualityConstraint(EigenSMatrix& Ai, double bi);
    void addInequalityConstraint(EigenSMatrix& Gi, double hi);

    void solve();

    void print() const;

    double getObjectiveValue() const;

    const EigenVector& getResult() const;

  private:

    int n_; // X is n by n matrix
    int m_; // Num Equality   Constraints
    int p_; // Num Inequality Constraints

    Parameters param_;

    int target_rank_;

    double rho_;

    double primal_infeasibility_;
    double primal_obj_;
    double runtime_;

    double b_l2_norm_;
    double b_inf_norm_;
    double C_inf_norm_;

    bool is_maximization_;

    /* ------------ CPU Data ------------ */
    EigenVector solution_;
    std::vector<double> b_;
    std::vector<double> h_;
    EigenSMatrix C_;
    std::vector<EigenSMatrix> A_;
    std::vector<EigenSMatrix> G_;

    /* ---------- CPU Function ---------- */
    void initParams();
    void printBanner() const;
    void printMetricRow() const;
    void printStartInfo(int target_rank, double initial_rho) const;
    void printFinish() const;
    void printProgress(
      int    alm_iter,
      int    lbfgs_iter,
      double rho,
      double primal_obj,
      double primal_inf,
      double runtime) const;
    /* ---------------------------------- */

    /* ------------ GPU Data ------------ */
    CudaFlattenMatrix<double>   d_C_; // const

    CudaFlattenMatrix<double>   d_UVT_workspace_;
    CudaFlattenMatrix<double>   d_R_;
    CudaFlattenMatrix<double>   d_Rnew_;

    CudaFlattenMatrix<double>   d_RGrad_;
    CudaFlattenMatrix<double>   d_RGradnew_;

    CudaFlattenMatrix<double>   d_LbfgsDirection_;
    CudaFlattenMatrix<double>   d_LbfgsWorkspace_;

    CudaFlattenMatrix<double>   d_X_;

    CudaVector<double>          d_b_; // const
    CudaVector<double>          d_h_; // const
   
    CudaVector<double>          d_linesearch_workspace_;

    CudaVector<double>          d_grad_workspace_vector_;
    CudaFlattenMatrix<double>   d_grad_workspace_matrix_;

    CudaSparseMatrix<double>    d_M_;
    CudaSparseMatrix<double>    d_MT_;

    std::vector<CudaFlattenMatrix<double>> d_A_;

    CudaVector<double>          d_lambda_;
    CudaVector<double>          d_ARRT_;
    CudaVector<double>          d_ARDT_;
    CudaVector<double>          d_ADDT_;

    CudaVector<double>          d_infeasibility_workspace_;

    std::shared_ptr<LbfgsNode>  lbfgs_head_;

    /* ---------- GPU Function ---------- */
    void solveALM();
    void initializeRho();
    void findTargetRank();
    void precomputeNorm();

    void precondition();

    void setInitialR(CudaFlattenMatrix<double>& R);

    void updateDual(
      const double rho,
      const CudaVector<double>& b,
      const CudaVector<double>& ARRT,
            CudaVector<double>& lambda);

    void computeUVT(
      const CudaFlattenMatrix<double>& U, 
      const CudaFlattenMatrix<double>& V, 
            CudaFlattenMatrix<double>& UVT);
    
    void computeConstraintValue(
      const CudaFlattenMatrix<double>& U,
      const CudaFlattenMatrix<double>& V,
            CudaVector<double>&        AUVT);

    double computeObjectiveValue(
      const CudaFlattenMatrix<double>& U,
      const CudaFlattenMatrix<double>& V);

    double lineSearchLbfgs(
      const double              rho,
      const double              CDDT,
      const double              CRDT,
      const CudaVector<double>& b,
      const CudaVector<double>& ARRT,
      const CudaVector<double>& ARDT,
      const CudaVector<double>& ADDT,
      const CudaVector<double>& lambda);

    void loadDataOnGPU();
    void makeSparseMatrixM();

    void computeGradient(
      const double                     rho,
      const CudaVector<double>&        b,
      const CudaVector<double>&        lambda,
      const CudaVector<double>&        ARRT,
      const CudaFlattenMatrix<double>& R,
      const CudaFlattenMatrix<double>& C,
            CudaFlattenMatrix<double>& Grad);

    int solveLbfgs(
      const double                     rho,
      const CudaVector<double>&        b,
      const CudaVector<double>&        lambda,
      const CudaFlattenMatrix<double>& C,
            CudaFlattenMatrix<double>& R);

    void computeLbfgsDirection(
      const int                        lbfgs_iter,
      const double                     rho,
      const CudaVector<double>&        b,
      const CudaVector<double>&        lambda,
      const CudaVector<double>&        ARRT,
      const CudaFlattenMatrix<double>& C,
      const CudaFlattenMatrix<double>& R,
      const CudaFlattenMatrix<double>& RGrad,
            CudaFlattenMatrix<double>& D); /* output */

    bool checkLbfgsConvergence(
      const double rho, 
      const CudaFlattenMatrix<double>& Grad);

    bool checkAlmConvergence(
      const CudaVector<double>& ARRT);

    /* ---------------------------------- */
};

}

#endif
