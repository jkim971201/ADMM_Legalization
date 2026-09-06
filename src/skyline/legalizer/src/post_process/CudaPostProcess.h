#ifndef LG_CUDA_POST_PROCESS_ENGINE_H
#define LG_CUDA_POST_PROCESS_ENGINE_H

#include "cuda_linalg/CudaVector.h"

namespace legalizer
{

class CudaDatabase;
class LGHyperParameters;

class CudaPostProcess
{
  public:

    CudaPostProcess(
      const std::shared_ptr<LGHyperParameters> hyperparameter,
            std::shared_ptr<CudaDatabase> database);

    void doPostProcess();

  private:

    int num_movable_cells_;
    int num_cut_rows_;

    void rotateCells();
    void stampCells();
    void makeCutRows();
    void insertCellsToCutRows();
    void sortCellsEachCutRow();
    void solveQP();

    cuda_linalg::CudaVector<int> grid_y_to_num_cut_rows_;
    cuda_linalg::CudaVector<int> grid_y_to_cut_row_offset_;

    cuda_linalg::CudaVector<int> cut_row_id_to_grid_y_;
    cuda_linalg::CudaVector<int> cut_row_id_to_grid_x_bgn_;
    cuda_linalg::CudaVector<int> cut_row_id_to_grid_x_end_;
    cuda_linalg::CudaVector<int> cut_row_id_to_cell_offset_;
    cuda_linalg::CudaVector<int> cut_row_id_to_num_cells_;

    cuda_linalg::CudaVector<int> pixel_id_to_cell_id_;
    cuda_linalg::CudaVector<int> pixel_id_to_cut_row_id_;

    cuda_linalg::CudaVector<int> cell_id_to_cut_row_id_;
    cuda_linalg::CudaVector<int> cell_grouped_by_cut_row_;

    std::shared_ptr<LGHyperParameters> params_;
    std::shared_ptr<CudaDatabase> database_;
};

}

#endif
