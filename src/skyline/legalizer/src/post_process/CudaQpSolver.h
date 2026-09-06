#ifndef CUDA_QP_SOLVER_H
#define CUDA_QP_SOLVER_H

#include "cuda_linalg/CudaVector.h"

namespace legalizer {
  class CudaDatabase;
  class LGHyperParameters;
}

namespace qp_solver
{

using namespace legalizer;

class CudaQpSolver
{
  public:

    CudaQpSolver(
      const std::shared_ptr<LGHyperParameters> hyperparameter,
            std::shared_ptr<CudaDatabase> database);

    void solve();

    void initializeSolver(
      const cuda_linalg::CudaVector<int>& cut_row_id_to_grid_x_bgn,
      const cuda_linalg::CudaVector<int>& cut_row_id_to_grid_x_end,
      const cuda_linalg::CudaVector<int>& cut_row_id_to_cell_offset,
      const cuda_linalg::CudaVector<int>& cut_row_id_to_num_cells,
      const cuda_linalg::CudaVector<int>& cell_grouped_by_cut_row);

  private:

    bool checkTermination(int iter);

    void updatePrimalX();
    void updatePrimalY();
    void updateDual();
    void computeConstraintValue();

    void siteAlignNaive();
    void siteAlignDP();

    float computePrimalObjective(const cuda_linalg::CudaVector<float>& primal) const;
    float computeDisplace(const cuda_linalg::CudaVector<float>& primal) const;

    float rho_;
    int num_movable_;
    int num_cut_row_;    // Num CutRows
    int num_constraint_; // Num Constraints (= NumCells + NumRows)

    std::shared_ptr<LGHyperParameters> params_;
    std::shared_ptr<CudaDatabase> database_;

    cuda_linalg::CudaVector<int>   global_id_to_lx_in_grid_;

    cuda_linalg::CudaVector<int>   local_id_to_global_id_;
    cuda_linalg::CudaVector<int>   local_id_to_pos_constraint_id_;
    cuda_linalg::CudaVector<int>   local_id_to_neg_constraint_id_;
    cuda_linalg::CudaVector<float> local_id_to_primal_; // cell lx_dbu in float
    cuda_linalg::CudaVector<float> local_id_to_original_lx_;

    cuda_linalg::CudaVector<float> constraint_id_to_residual_;
    cuda_linalg::CudaVector<float> constraint_id_to_primal_product_;
    cuda_linalg::CudaVector<float> constraint_id_to_slack_;
    cuda_linalg::CudaVector<float> constraint_id_to_dual_;
    cuda_linalg::CudaVector<float> constraint_id_to_bound_;
    cuda_linalg::CudaVector<int>   constraint_id_to_pos_local_cell_id_;
    cuda_linalg::CudaVector<int>   constraint_id_to_neg_local_cell_id_;

    cuda_linalg::CudaVector<int>   cut_row_id_to_num_cells_;
    cuda_linalg::CudaVector<int>   cut_row_id_to_cell_offset_;
    cuda_linalg::CudaVector<int>   cut_row_id_to_num_constraint_;
    cuda_linalg::CudaVector<int>   cut_row_id_to_constraint_offset_;
    cuda_linalg::CudaVector<int>   cut_row_id_to_grid_x_bgn_;
    cuda_linalg::CudaVector<int>   cut_row_id_to_grid_x_end_;

    cuda_linalg::CudaVector<float> thomas_precomputed_c_prime_;
    cuda_linalg::CudaVector<float> thomas_workspace_c_prime_;
    cuda_linalg::CudaVector<float> thomas_workspace_d_prime_;
    cuda_linalg::CudaVector<float> thomas_workspace_d_in_;
};

}

#endif
