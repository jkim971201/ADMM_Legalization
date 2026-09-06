#ifndef LG_CUDA_ADMM_LEGALIZER_H
#define LG_CUDA_ADMM_LEGALIZER_H

#include <map>
#include <deque>
#include "cuda_linalg/CudaVector.h"

namespace legalizer {
  class Grid;
  class CudaDatabase;
  class LGHyperParameters;
}

namespace admm_legalizer
{

using namespace legalizer;

class CudaAdmmLegalizer
{
  public:

    CudaAdmmLegalizer(
      const std::vector<int>& direction_default,
      const std::vector<int>& direction_horizontal,
      const std::vector<int>& direction_vertical,
      const std::vector<int>& direction_wide,
      const std::map<int, int>& height_distribution,
      const std::shared_ptr<Grid> grid,
      const std::shared_ptr<LGHyperParameters> hyperparameter,
            std::shared_ptr<CudaDatabase> database);

    void solve();

  private:

    void importFromGlobalDatabase();
    void exportToGlobalDatabase();

    void applyAdmm();
      void updatePartition(int iter);
        void paintPixels(int offset_x, int offset_y);

      std::pair<int, float> examineIteration();

      float computeInverseAugmentedDemand();
      void initializeRho();
      void initializeDual();
      void doPrimalUpdate(int iter);
      void doDualUpdate();
      void updateDualStepSize();

      void printProgress(int iter, int ovf, float disp) const;

      bool detectSlowConvergence(int iter, int ovf, int min_ovf_so_far);

      void doPrimalUpdateForSlowConvergence(int iter);
      void doFinalRefine();

      void executeSorting(
        const cuda_linalg::CudaVector<int>& lx_in_grid,
        const cuda_linalg::CudaVector<int>& ly_in_grid,
        const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
              cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition);

      void sortEachPartition(
        const int num_threads_per_block,
        const cuda_linalg::CudaVector<int>& lx_in_grid,
        const cuda_linalg::CudaVector<int>& ly_in_grid,
        const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
              cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition);

      void sortEachPartitionWithMergeSort(
        const int num_threads_per_block,
        const cuda_linalg::CudaVector<int>& lx_in_grid,
        const cuda_linalg::CudaVector<int>& ly_in_grid,
        const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
              cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition);

      void doBitonicSort(
        const int num_threads_per_block,
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
              cuda_linalg::CudaVector<int>& cell_ovf_grouped_by_partition,  // Key
              cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition);  // Value

      void doMergeSort(
        const int num_threads_per_block,
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
              cuda_linalg::CudaVector<int>& cell_ovf_grouped_by_partition,  // Key
              cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition);  // Value

      void doCubSort(
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
              cuda_linalg::CudaVector<int>& cell_ovf_grouped_by_partition,  // Key
              cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition);  // Value

      void computeKeyArrayForSorting(
        const cuda_linalg::CudaVector<int>& lx_in_grid,
        const cuda_linalg::CudaVector<int>& ly_in_grid,
        const cuda_linalg::CudaVector<int>& partition_id_to_num_cell,
        const cuda_linalg::CudaVector<int>& partition_id_to_offset,
        const cuda_linalg::CudaVector<int>& cell_id_grouped_by_partition,
              cuda_linalg::CudaVector<int>& cell_key_grouped_by_partition);

    void writePartition(int iter) const;

    void writeRecord() const;

    int num_cells_;
    int site_width_;
    int row_height_;

    int x_size_;
    int y_size_;

    int num_partition_col_;
    int num_partition_row_;

    int max_num_cell_in_partition_;

    int core_lx_;
    int core_ly_;

    int num_directions_padded_;

    int min_ovf_so_far_;
    int last_iter_slow_convergence_detected_;

    // Metrics
    std::deque<int> ovf_dequeue_;
    std::vector<std::pair<int, float>> iter_record_;

    // ADMM parameters 
    float rho_;
    float dual_step_size_;
    float initial_dual_;

    std::shared_ptr<LGHyperParameters> params_;
    std::shared_ptr<CudaDatabase> database_;

    cuda_linalg::CudaVector<int>   cell_random_priority_;

    cuda_linalg::CudaVector<int>   direction_default_;
    cuda_linalg::CudaVector<int>   direction_horizontal_;
    cuda_linalg::CudaVector<int>   direction_vertical_;
    cuda_linalg::CudaVector<int>   direction_wide_;

    cuda_linalg::CudaVector<int>   sequence_partition_;
    cuda_linalg::CudaVector<int>   partition_id_to_offset_;
    cuda_linalg::CudaVector<int>   partition_id_to_num_cells_;
    cuda_linalg::CudaVector<int>   cell_id_to_partition_id_;
    cuda_linalg::CudaVector<int>   cell_id_grouped_by_partition_;
    cuda_linalg::CudaVector<int>   cell_key_grouped_by_partition_;

    cuda_linalg::CudaVector<int>   cell_id_to_is_ovf_;
    cuda_linalg::CudaVector<int>   ovf_cells_;
    cuda_linalg::CudaVector<int>   cell_compaction_workspace_;

    cuda_linalg::CudaVector<int>   k_height_to_num_cells_;

    cuda_linalg::CudaVector<int>   cell_lx_in_grid_;
    cuda_linalg::CudaVector<int>   cell_ly_in_grid_;

    cuda_linalg::CudaVector<float> cell_demand_plus_area_;
    cuda_linalg::CudaVector<float> cell_inverse_augmented_demand_;

    cuda_linalg::CudaVector<float> pixel_dual_;
    cuda_linalg::CudaVector<int>   pixel_color_; // color = partition_id
};

}

#endif
