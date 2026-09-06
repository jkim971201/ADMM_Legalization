#ifndef LG_CUDA_DATABASE_H
#define LG_CUDA_DATABASE_H

#include "cuda_linalg/CudaVector.h"

namespace legalizer
{

class Grid;
class LGRegion;
class LGCell;
class LGCellTechInfo;
class LGHyperParameters;
struct LGVioStat;

class CudaDatabase
{
  public:

    CudaDatabase(
      const std::vector<std::shared_ptr<LGRegion>>& regions,
      const std::shared_ptr<LGHyperParameters> param, 
      const std::shared_ptr<Grid> grid);

    LGVioStat checkPlace();

    void copyEdgeSpacingRules(const std::vector<int>& edge_spacing_rules);

    void copyFromHostCells(
      const std::vector<LGCell>& cells,
      const std::vector<std::shared_ptr<LGCellTechInfo>>& tech_infos);

    void copyToHostCells(std::vector<LGCell>& cells);

    void commitPlaceInGrid(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid);

    void updateUsageInDbu(
      const cuda_linalg::CudaVector<int>& lx_in_dbu,
      const cuda_linalg::CudaVector<int>& ly_in_dbu);

    void updateUsageInGrid(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid);

    void updatePixelEdgeMapInGrid(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid);

    int getCoreLx() const { return core_lx_; }
    int getCoreLy() const { return core_ly_; }
    int getXSize() const { return x_size_; }
    int getYSize() const { return y_size_; }
    int getNumCells() const { return num_cells_; }
    int getNumPixels() const { return x_size_ * y_size_; }

    int getSiteWidth() const { return site_width_; }
    int getRowHeight() const { return row_height_; }

    int computeOverflow() const;

    float computeDisplace() const;
    float computeDisplace(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid) const;

    const cuda_linalg::CudaVector<int>& getPinBgn() const { return pin_id_to_bgn_; }
    const cuda_linalg::CudaVector<int>& getPinEnd() const { return pin_id_to_end_; }
    const cuda_linalg::CudaVector<int>& getPinLayer() const { return pin_id_to_layer_; }
    const cuda_linalg::CudaVector<int>& getPinOffsets() const { return cell_macro_id_to_pin_offset_; }
    const cuda_linalg::CudaVector<int>& getCellMacroID() const { return cell_id_to_cell_macro_id_; }

    const cuda_linalg::CudaVector<int>& getCellWidthInGrid() const { return cell_id_to_width_in_grid_; }
    const cuda_linalg::CudaVector<int>& getCellHeightInGrid() const { return cell_id_to_height_in_grid_; }

    const cuda_linalg::CudaVector<int>& getCellHasGroundAtBottom() const { return cell_id_to_has_ground_at_bottom_; }

    const cuda_linalg::CudaVector<int>& getCellLEdgeType() const { return cell_id_to_ledge_type_; }
    const cuda_linalg::CudaVector<int>& getCellREdgeType() const { return cell_id_to_redge_type_; }
    const cuda_linalg::CudaVector<int>& getEdgeSpacingRules() const { return edge_spacing_rules_; }

    const cuda_linalg::CudaVector<int>& getCellGroupId() const { return cell_id_to_group_index_; }

    const cuda_linalg::CudaVector<int>& getCellOriginalLxInDbu() const { return cell_id_to_original_lx_in_dbu_; }
    const cuda_linalg::CudaVector<int>& getCellOriginalLyInDbu() const { return cell_id_to_original_ly_in_dbu_; }

    const cuda_linalg::CudaVector<int>& getCellLxInGrid() const { return cell_id_to_lx_in_grid_; }
    const cuda_linalg::CudaVector<int>& getCellLyInGrid() const { return cell_id_to_ly_in_grid_; }

    const cuda_linalg::CudaVector<int>& getCellNeedRotation() const { return cell_id_to_need_rotation_; }
          cuda_linalg::CudaVector<int>& getCellNeedRotation()       { return cell_id_to_need_rotation_; }

    const cuda_linalg::CudaVector<int>& getCellNeedFlip() const { return cell_id_to_need_flip_; }
          cuda_linalg::CudaVector<int>& getCellNeedFlip()       { return cell_id_to_need_flip_; }

    const cuda_linalg::CudaVector<int>& getGridYToIsVddUp() const { return grid_y_to_is_vdd_up_; }

    const cuda_linalg::CudaVector<int>& getPixelValid() const { return pixel_id_to_valid_; }

    const cuda_linalg::CudaVector<int>& getPixelGroupId() const { return pixel_id_to_group_id_; }

    const cuda_linalg::CudaVector<int>& getPixelRegionId() const { return pixel_id_to_region_id_; }

    const cuda_linalg::CudaVector<int>& getPixelPowerMetalLayer() const { return pixel_id_to_power_metal_layer_; }

    const cuda_linalg::CudaVector<int>& getPixelUsage() const { return pixel_id_to_usage_; }
          cuda_linalg::CudaVector<int>& getPixelUsage()       { return pixel_id_to_usage_; }

    const cuda_linalg::CudaVector<int>& getPixelShape() const { return pixel_id_to_shape_; }

          cuda_linalg::CudaVector<int>& getPixelNumLEdgeType1() { return pixel_id_to_num_ledge_type1_; }
          cuda_linalg::CudaVector<int>& getPixelNumLEdgeType2() { return pixel_id_to_num_ledge_type2_; }
          cuda_linalg::CudaVector<int>& getPixelNumREdgeType1() { return pixel_id_to_num_redge_type1_; }
          cuda_linalg::CudaVector<int>& getPixelNumREdgeType2() { return pixel_id_to_num_redge_type2_; }

    const cuda_linalg::CudaVector<int>& getRegionGridXMin() const { return region_id_to_grid_x_min_; }
    const cuda_linalg::CudaVector<int>& getRegionGridYMin() const { return region_id_to_grid_y_min_; }
    const cuda_linalg::CudaVector<int>& getRegionGridXMax() const { return region_id_to_grid_x_max_; }
    const cuda_linalg::CudaVector<int>& getRegionGridYMax() const { return region_id_to_grid_y_max_; }

  private:

    void writeViolationMap(const cuda_linalg::CudaVector<int>& access_vio) const;

    void copyGridInfo(const std::vector<std::shared_ptr<LGRegion>>& regions);

    void markPinLayerToGrid(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid,
            cuda_linalg::CudaVector<int>& pixel_id_to_m2_pin);

    void convertDbuToGrid(
      const cuda_linalg::CudaVector<int>& lx_in_dbu,
      const cuda_linalg::CudaVector<int>& ly_in_dbu,
            cuda_linalg::CudaVector<int>& lx_in_grid,
            cuda_linalg::CudaVector<int>& ly_in_grid);

    void convertGridToDbu(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid,
            cuda_linalg::CudaVector<int>& lx_in_dbu,
            cuda_linalg::CudaVector<int>& ly_in_dbu);

    void countPinVio(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid,
            cuda_linalg::CudaVector<int>& cell_id_to_num_pin_short_vio,
            cuda_linalg::CudaVector<int>& cell_id_to_num_pin_access_vio);

    void countEdgeSpacingVio(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid,
            cuda_linalg::CudaVector<int>& cell_id_to_num_edge_spacing_vio);

    void countRowOrientVio(
      const cuda_linalg::CudaVector<int>& lx_in_grid,
      const cuda_linalg::CudaVector<int>& ly_in_grid,
            cuda_linalg::CudaVector<int>& cell_id_to_num_row_orient_vio);

    void diagnosePixels();

    int num_cells_;
    int x_size_;
    int y_size_;

    int site_width_;
    int row_height_;

    int core_lx_;
    int core_ly_;

    std::shared_ptr<LGHyperParameters> param_;
    std::shared_ptr<Grid> grid_;

    cuda_linalg::CudaVector<int> pixel_id_to_usage_;
    cuda_linalg::CudaVector<int> pixel_id_to_valid_;
    cuda_linalg::CudaVector<int> pixel_id_to_group_id_;

    cuda_linalg::CudaVector<int> pixel_id_to_shape_;

    cuda_linalg::CudaVector<int> pixel_id_to_num_ledge_type1_;
    cuda_linalg::CudaVector<int> pixel_id_to_num_ledge_type2_;
    cuda_linalg::CudaVector<int> pixel_id_to_num_redge_type1_;
    cuda_linalg::CudaVector<int> pixel_id_to_num_redge_type2_;

    cuda_linalg::CudaVector<int> pixel_id_to_power_metal_layer_;;

    cuda_linalg::CudaVector<int> pixel_id_to_region_id_;

    cuda_linalg::CudaVector<int> region_id_to_grid_x_min_;
    cuda_linalg::CudaVector<int> region_id_to_grid_y_min_;
    cuda_linalg::CudaVector<int> region_id_to_grid_x_max_;
    cuda_linalg::CudaVector<int> region_id_to_grid_y_max_;

    cuda_linalg::CudaVector<int> grid_y_to_is_vdd_up_;

    cuda_linalg::CudaVector<int> cell_id_to_has_ground_at_bottom_;

    cuda_linalg::CudaVector<int> edge_spacing_rules_;
    cuda_linalg::CudaVector<int> cell_id_to_ledge_type_;
    cuda_linalg::CudaVector<int> cell_id_to_redge_type_;

    cuda_linalg::CudaVector<int> cell_id_to_cell_macro_id_;

    cuda_linalg::CudaVector<int> cell_id_to_num_short_vio_;
    cuda_linalg::CudaVector<int> cell_id_to_num_access_vio_;
    cuda_linalg::CudaVector<int> cell_id_to_num_edge_spacing_vio_;

    cuda_linalg::CudaVector<int> cell_macro_id_to_pin_offset_;
    cuda_linalg::CudaVector<int> pin_id_to_bgn_;
    cuda_linalg::CudaVector<int> pin_id_to_end_;
    cuda_linalg::CudaVector<int> pin_id_to_layer_;

    cuda_linalg::CudaVector<int> cell_id_to_original_lx_in_dbu_;
    cuda_linalg::CudaVector<int> cell_id_to_original_ly_in_dbu_;

    cuda_linalg::CudaVector<int> cell_id_to_lx_in_grid_;
    cuda_linalg::CudaVector<int> cell_id_to_ly_in_grid_;

    cuda_linalg::CudaVector<int> cell_id_to_lx_in_dbu_;
    cuda_linalg::CudaVector<int> cell_id_to_ly_in_dbu_;

    cuda_linalg::CudaVector<int> cell_id_to_width_in_grid_;
    cuda_linalg::CudaVector<int> cell_id_to_height_in_grid_;

    cuda_linalg::CudaVector<int> cell_id_to_group_index_;

    cuda_linalg::CudaVector<int> cell_id_to_need_rotation_;
    cuda_linalg::CudaVector<int> cell_id_to_need_flip_;
};

}

#endif
