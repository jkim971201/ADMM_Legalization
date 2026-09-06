#include <fstream>
#include <ranges>
#include <cassert>

#include "legalizer/Legalizer.h"
#include "objects/LGCell.h"
#include "objects/LGCellTechInfo.h"
#include "objects/LGRegion.h"

#include "db/dbTypes.h"
#include "db/dbDatabase.h"
#include "db/dbDesign.h"
#include "db/dbInst.h"
#include "db/dbMacro.h"
#include "db/dbTech.h"
#include "db/dbLayer.h"
#include "db/dbNet.h"
#include "db/dbRow.h"
#include "db/dbDie.h"
#include "db/dbRect.h"
#include "db/dbRegion.h"
#include "db/dbBlockage.h"

#include "util/Chrono.h"

#include "Grid.h"

namespace legalizer
{

void
Legalizer::importDB()
{
  auto import_start = util::getChronoNow();

  // Step 1. Import Chip Boundary
  const auto design = dbDatabase_->getDesign();
  core_lx_ = design->coreLx();
  core_ly_ = design->coreLy();
  core_ux_ = design->coreUx();
  core_uy_ = design->coreUy();

  // Step 2. Import Cells 
  importCells();

  // Step 3. Create Grid
  createGrid();

  // Step 4. Import Region/Group Info
  for(int group_index = 0; const auto& db_group_ptr : design->getGroups())
    db_group_to_index_[db_group_ptr] = group_index++;

  importBlockageRegions();
  importGroupRegions();
  importNonGroupRegions();
  assignGroupMembers();

  findRegionAdjacency();

  // Step 5. Mark Region
  markRegion();
  handleRegionSegments();

  // Step 6. Import Technology Info
  importTechInfo();

  // Step 7. Compute Original HPWL
  original_hpwl_ = computeHpwlFromDB();

  // Read placement.constraints file of ICCAD2017 benchmarks.
  if(!constraint_.empty())
    parseICCAD17Constraint();

  printDesignInfo();

  printf("importDB           finished (takes %5.2f s)\n", util::evalTime(import_start));
}

void
Legalizer::importCells()
{
  const auto design                    = dbDatabase_->getDesign();
  const std::vector<dbInst*>& db_insts = design->getInsts();
  const std::vector<dbRow*>& db_rows   = design->getRows();
  const int site_width                 = db_rows.back()->siteWidth(); // TODO: Fix this.
  const int row_height                 = db_rows.back()->dy();
  // <- Assume uniform row config

  std::unordered_map<std::string, std::pair<int ,int>> ispd2015_size;
  bool modifiedISPD2015 = (size_file_.empty() == true) ? false : true;
  if(modifiedISPD2015 == true)
  {
    printf("Running with Modifed ISPD 2015 Benchmark\n");
    parseSizeFile(ispd2015_size);
  }

  sum_cell_area_  = 0;
  sum_fixed_area_ = 0;
  sum_movable_area_ = 0;
  for(auto dbinst_ptr : db_insts)
  {
    if(dbinst_ptr->dx() == 0 || dbinst_ptr->dy() < row_height)
      continue;

    int64_t area_this_cell = int64_t(dbinst_ptr->area());
    sum_cell_area_ += area_this_cell;
    if(dbinst_ptr->isMacro())
    {
      sum_fixed_area_ += area_this_cell;
      if(!dbinst_ptr->isFixed())
      {
        printf("Cannot handle movable macro yet...\n");
        exit(1);
      }
      else
      {
        LGRect rect_of_fixed_macro(dbinst_ptr);
        fixed_macro_rects_.push_back(rect_of_fixed_macro);
      }
    }
    else
    {
      if(!dbinst_ptr->isFixed())
      {
        sum_movable_area_ += area_this_cell;
        LGCell new_cell(dbinst_ptr);

        if(modifiedISPD2015 == true)
        {
          auto size_itr = ispd2015_size.find(dbinst_ptr->name());
          if(size_itr != ispd2015_size.end())
          {
            const auto [new_width, new_height] = size_itr->second;
            new_cell.setWidth(new_width * site_width);
            new_cell.setHeight(new_height * row_height);
          }
        }

        cells_.push_back(new_cell);
        height_distribution_[new_cell.getHeight()] += 1;
      }
      else
      {
        sum_fixed_area_ += area_this_cell;
        LGRect rect_of_fixed_stdcell(dbinst_ptr);
        fixed_stdcell_rects_.push_back(rect_of_fixed_stdcell);
      }
    }
  }

  // Make dbInst2LGCell hash map
  for(auto& cell : cells_)
    db_inst_to_lg_cell_[cell.getDbInst()] = &cell;

//  for(int cell_id = 0; cell_id < cells_.size(); cell_id++)
//  {
//    const auto& cell = cells_[cell_id];
//    if(cell.getName() == "o1283778")
//      printf("CellName: %s -> CellID: %d\n", cell.getName().data(), cell_id);
//    if(cell.getName() == "o1283845")
//      printf("CellName: %s -> CellID: %d\n", cell.getName().data(), cell_id);
//    if(cell.getName() == "o1284141")
//      printf("CellName: %s -> CellID: %d\n", cell.getName().data(), cell_id);
//    if(cell.getName() == "o1284206")
//      printf("CellName: %s -> CellID: %d\n", cell.getName().data(), cell_id);
//    if(cell.getName() == "o1284365")
//      printf("CellName: %s -> CellID: %d\n", cell.getName().data(), cell_id);
//  }
}

void
Legalizer::assignGroupMembers()
{
  const auto& groups = dbDatabase_->getDesign()->getGroups();
  for(const auto& db_group_ptr : groups)
  {
    const int group_index = db_group_to_index_[db_group_ptr];
    const auto& members_this_group = db_group_ptr->getMembers();
    for(const auto db_inst : members_this_group)
    {
      LGCell* cell = db_inst_to_lg_cell_[db_inst];
      cell->setGroupIndex(group_index);
    }
  }
}

void
Legalizer::createGrid()
{
  grid_           = std::make_shared<Grid>(dbDatabase_);
  site_width_     = grid_->getSiteWidth();
  row_height_     = grid_->getRowHeight();
  num_total_rows_ = grid_->getYSize(); // only count for unique rows (not cut rows)
}

void
Legalizer::markRegion()
{
  for(const auto region : regions_)
  {
    if(region->getType() == RegionType::BLOCKAGE)
      grid_->markBlockageRegion(region);
    else if(region->getType() == RegionType::FENCE)
      grid_->markGroupRegion(region);
    else if(region->getType() == RegionType::DEFAULT)
      grid_->markNonGroupRegion(region);
  }
}

void
Legalizer::importBlockageRegions()
{
  for(const auto& rect : fixed_macro_rects_)
  {
    // There can be blocks out of core bbox
    int region_lx = std::max(core_lx_, rect.getLx());
    int region_ly = std::max(core_ly_, rect.getLy());
    int region_ux = std::min(core_ux_, rect.getUx());
    int region_uy = std::min(core_uy_, rect.getUy());

    auto new_region 
      = std::make_shared<LGRegion>(regions_.size(), // region_id
          region_lx, region_ly, region_ux, region_uy,
          k_blockage_group_index, RegionType::BLOCKAGE);
    addNewRegion(new_region);
  }

  for(const auto& rect : fixed_stdcell_rects_)
  {
    // There can be blocks out of core bbox
    int region_lx = std::max(core_lx_, rect.getLx());
    int region_ly = std::max(core_ly_, rect.getLy());
    int region_ux = std::min(core_ux_, rect.getUx());
    int region_uy = std::min(core_uy_, rect.getUy());

    auto new_region 
      = std::make_shared<LGRegion>(regions_.size(), // region_id
          region_lx, region_ly, region_ux, region_uy,
          k_blockage_group_index, RegionType::BLOCKAGE);
    addNewRegion(new_region);
  }
}

void
Legalizer::importGroupRegions()
{
  const auto design = dbDatabase_->getDesign();
  const auto db_regions = design->getRegions();

  for(auto db_region_ptr : db_regions)
  {
    dbGroup* db_group_ptr = db_region_ptr->getGroup();
    const int group_index = db_group_to_index_[db_group_ptr];
    for(auto& box : db_region_ptr->getRegionBoxes())
    {
      // There can be group regions out of core bbox
      int region_lx = std::max(core_lx_, box.lx());
      int region_ly = std::max(core_ly_, box.ly());
      int region_ux = std::min(core_ux_, box.ux());
      int region_uy = std::min(core_uy_, box.uy());

      auto new_region 
        = std::make_shared<LGRegion>(regions_.size(), // region_id
            region_lx, region_ly, region_ux, region_uy,
            group_index, RegionType::FENCE);
      addNewRegion(new_region);
    }
  }
}

void
Legalizer::importNonGroupRegions()
{
  // For linesweep algorithm
  struct Event 
  { 
    int rect_id, y;
    bool is_enter;
  };

  std::vector<Event> events;

  std::vector<const LGRect*> rects;
  for(auto region : regions_)
  {
    const RegionType type = region->getType();
    if(type == RegionType::FENCE or type == RegionType::BLOCKAGE)
      rects.push_back(region.get());
  }

  for(int rect_id = 0; auto rect : rects) // <- C++ 20 feature?
  {
    events.push_back({rect_id, rect->getLy(), true});  
    events.push_back({rect_id, rect->getUy(), false}); 
    rect_id++;
  }

  std::sort(events.begin(), events.end(), [] (const Event& a, const Event& b) 
    {
      if(a.y == b.y) return !a.is_enter && b.is_enter;
      else return a.y < b.y;
    });

  std::set<int> active_rects;

  // we need "lx <-> regions table" to merge regions that have same width and lx
  std::unordered_map<int, std::vector<std::shared_ptr<LGRegion>>> lx_to_regions;

  auto create_region = [&] (int lx, int ly, int ux, int uy)
  {
    // if region size is too small, it's not valid.
    if((ux - lx) < site_width_ or (uy - ly) < row_height_)
      return;
    else
    {
      const int width = ux - lx;
      auto& regions_this_lx = lx_to_regions[lx];

      bool merge_done = false;
      for(auto& region : regions_this_lx)
      {
        // If having same width and adjacent (which means its uy equals to new region's ly)
        if(region->getWidth() == width && region->getUy() == ly)
        {
          int new_height = uy - region->getLy();
          region->setHeight(new_height);
          merge_done = true;
          break;
        }
      }

      if(merge_done == false)
      {
        //printf("LxLy (%7d, %7d) UxUy (%7d, %7d)\n", lx, ly, ux, uy);
        auto new_region = std::make_shared<LGRegion>(regions_.size(), // region_id
          lx, ly, ux, uy, k_default_group_index, RegionType::DEFAULT);
        regions_this_lx.push_back(new_region);
        addNewRegion(new_region);
      }
    }
  };

  auto make_partial_nongroup_regions = [&] (int ly, int uy)
  {
    std::vector<std::pair<int, int>> intervals;
    for(int rect_id : active_rects)
    {
      auto rect = rects[rect_id];
      intervals.push_back({rect->getLx(), rect->getUx()});
    }

    std::sort(intervals.begin(), intervals.end());
    // Default behavior -> ascending order of the first element
    // (use second element for tie-breaking)
    
    int lx = core_lx_;
    for(const auto& [interval_x1, interval_x2] : intervals)
    {
      create_region(lx, ly, interval_x1, uy);
      lx = interval_x2;
    }

    // rightmost partial region
    create_region(lx, ly, core_ux_, uy);
  };

  int ly = core_ly_;
  for(const auto& e : events) 
  {
    make_partial_nongroup_regions(ly, e.y);
    if(e.is_enter == true) 
      active_rects.insert(e.rect_id);
    else 
      active_rects.erase(e.rect_id);
    ly = e.y;
  }

  create_region(core_lx_, ly, core_ux_, core_uy_);
}

void
Legalizer::importTechInfo()
{
  const auto db_tech = dbDatabase_->getTech();
  const auto db_design = dbDatabase_->getDesign();

  importLayers(db_tech->getLayers());

  importEdgeSpcingRules(db_tech->getDbu(), db_tech->getEdgeSpacingRules());

  importCellMacro(db_tech->getMacros());

  importPowerMetals(db_design->getSNets());
}

void
Legalizer::importLayers(const std::vector<dbLayer*>& db_layers)
{
  // We will assign index only for routing_layers
  // (exclusing VIA or other cut/slice layers)
  std::vector<const dbLayer*> routing_layers;
  for(const auto& db_layer : db_layers)
  {
    if(db_layer->type() == RoutingType::ROUTING)
      routing_layers.push_back(db_layer);
  }

  std::sort(routing_layers.begin(), routing_layers.end(),
    [] (auto l1, auto l2) { return l1->index() < l2->index(); });

  for(int layer_index = 0; layer_index < routing_layers.size(); layer_index++)
    db_layer_to_index_[routing_layers[layer_index]] = layer_index;
}

void
Legalizer::importEdgeSpcingRules(int dbu, const std::vector<std::tuple<int, int, double>>& edge_spacing_rules)
{
  // Import Edge Spacing Rule into 1-D array (vector)
  // 0 : Empty
  // 1 : Empty
  // 1 + 1 : Spacing between type1 and type1
  // 1 + 2 : Spacing between type1 and type2
  // 2 + 2 : Spacing between type2 and type2
  // 5:  Empty
  // <- I know this is stupid, but I cannot think of better idea for array-based representation.
  // (since CUDA kernel only accepts array-style input)
  int max_type = 0;
  for(const auto& [type1, type2, spacing] : edge_spacing_rules)
  {
    max_type = std::max(max_type, type1);
    max_type = std::max(max_type, type2);
  }

  edge_spacing_rules_.resize(max_type * 2 + 1, 0);
  for(const auto& [type1, type2, spacing] : edge_spacing_rules)
  {
    int spacing_in_dbu = static_cast<int>(spacing * static_cast<float>(dbu));
    int spacing_in_grid = spacing_in_dbu / site_width_;
    edge_spacing_rules_[type1 + type2] = spacing_in_grid;
  }
}

void
Legalizer::importCellMacro(const std::vector<dbMacro*>& db_macros)
{
  // Create CellTechInfo
  int macro_id = 0;
  for(const auto& db_macro : db_macros)
  {
    // skip for macro block
    if(db_macro->macroClass() == MacroClass::BLOCK)
      continue;

    //printf("dbMacro: %s ", db_macro->name().c_str());
    //printf("Size (%d, %d)\n", db_macro->sizeX(), db_macro->sizeY());
    std::shared_ptr<LGCellTechInfo> new_tech_info = std::make_shared<LGCellTechInfo>(macro_id, db_macro);
    auto& pin_shapes = new_tech_info->getPinShapes();

    int cell_w_in_grid = db_macro->sizeX() / site_width_;
    int cell_h_in_grid = db_macro->sizeY() / row_height_;

    // We import pin shape here because we need site_width and row_height 
    // to flatten pin coordinates using sw and rh.
    const auto db_mterms = db_macro->getMTerms();
    for(const auto& mterm : db_mterms)
    {
      // Ignore POWER / GROUND pin
      if(mterm->usage() == PinUsage::POWER or mterm->usage() == PinUsage::GROUND)
        continue;

      const auto& ports_this_mterm = mterm->ports();

      int bbox_lx = std::numeric_limits<int>::max();
      int bbox_ly = std::numeric_limits<int>::max();
      int bbox_ux = 0;
      int bbox_uy = 0;
  
      int layer_index = -1;
      for(const auto& port : ports_this_mterm)
      {
        const std::vector<std::pair<int, int>>& shape = port->getShape();
        for(const auto& [offsetX_dbu, offsetY_dbu] : shape)
        {
          bbox_lx = std::min(bbox_lx, offsetX_dbu);
          bbox_ly = std::min(bbox_ly, offsetY_dbu);
          bbox_ux = std::max(bbox_ux, offsetX_dbu);
          bbox_uy = std::max(bbox_uy, offsetY_dbu);
        }
  
        // NOTE: we assume every port has same layer index
        layer_index = db_layer_to_index_[port->layer()];
      }

      if(layer_index == -1)
        continue;

      int bbox_lx_in_grid = bbox_lx / site_width_;
      int bbox_ly_in_grid = bbox_ly / row_height_;
      int bbox_ux_in_grid = std::min(cell_w_in_grid - 1, bbox_ux / site_width_);
      int bbox_uy_in_grid = std::min(cell_h_in_grid - 1, bbox_uy / row_height_);

      assert(bbox_lx_in_grid >= 0);
      assert(bbox_ly_in_grid >= 0);
      assert(bbox_ux_in_grid < cell_w_in_grid);
      assert(bbox_uy_in_grid < cell_h_in_grid);
  
      int bgn = bbox_lx_in_grid * cell_h_in_grid + bbox_ly_in_grid;
      int end = bbox_ux_in_grid * cell_h_in_grid + bbox_uy_in_grid;
  
      pin_shapes.push_back({layer_index, bgn, end});
      //printf("Layer: %d Bgn: %d End: %d\n", layer_index, bgn, end);
    }

    cell_tech_infos_.push_back(new_tech_info);
    db_macro_to_tech_info_[db_macro] = new_tech_info;
    macro_id++;
  }

  // Link CellTechInfo for each LGCells
  //int cell_id = 0;
  for(auto& cell : cells_)
  {
    auto db_macro = cell.getDbInst()->macro();
    auto& cell_tech_info = db_macro_to_tech_info_[db_macro];
    cell.setTechInfo(cell_tech_info);
    //if(db_macro->name() == "ms00f80" or db_macro->name() == "oa22f01")
    //  printf("Cell: %s ID: %d MacroID: %d\n", 
    //    cell.getName().data(), cell_id, cell_tech_info->getID());
    //cell_id++;
  }
}

void
Legalizer::importPowerMetals(const std::vector<dbNet*>& special_nets)
{
  // Mark Power Metals
  for(const auto& snet : special_nets)
  {
    const auto& segments = snet->getWire()->getSegments();
    for(const auto& seg : segments)
    {
      const dbLayer* seg_layer = seg->layer();
      if(seg_layer->type() != RoutingType::ROUTING)
        continue;

      int layer_index = db_layer_to_index_[seg_layer];
      if(layer_index == 0 or layer_index > 2) // Ignore M1 Metal and M3 Metal
        continue;

      int x1 = seg->startX();
      int x2 = seg->endX();
      if(x1 > x2) 
        std::swap(x1, x2);

      int y1 = seg->startY();
      int y2 = seg->endY();
      if(y1 > y2) 
        std::swap(y1, y2);

      // This is via?
      if(x1 == x2 and y1 == y2)
        continue;

      // we apply X2 to avoid spacing violation
      // TODO: Fix this hard code
      int width = seg->width();

      int seg_lx, seg_ly, seg_ux, seg_uy;
      seg_lx = std::max(core_lx_, x1 - width / 2);
      seg_ux = std::min(core_ux_, x2 + width / 2);
      seg_ly = std::max(core_ly_, y1 - width / 2);
      seg_uy = std::min(core_uy_, y2 + width / 2);

      //printf("Seg (%d, %d) - (%d, %d)\n", seg_lx, seg_ly, seg_ux, seg_uy);
      LGRect metal_rect(seg_lx, seg_ly, seg_ux - seg_lx, seg_uy - seg_ly);
      grid_->markPowerMetalLayer(layer_index, &metal_rect);
    }
  }
}

void
Legalizer::exportDB() const
{
  auto export_start = util::getChronoNow();

  for(const auto& cell : cells_)
  {
    int new_lx = cell.getLx();
    int new_ly = cell.getLy();
    auto db_inst_ptr = cell.getDbInst();
    db_inst_ptr->setLocation(new_lx, new_ly);
    db_inst_ptr->setOrient(cell.getOrient());
  }

  printf("exportDB           finished (takes %5.2f s)\n", util::evalTime(export_start));
}

void
Legalizer::findRegionAdjacency()
{
  const int num_regions = static_cast<int>(regions_.size());
  for(int i = 0; i < num_regions - 1; i++)
  {
    auto& region_i = regions_[i];
    const int lx_i = region_i->getLx();
    const int ly_i = region_i->getLy();
    const int ux_i = region_i->getUx();
    const int uy_i = region_i->getUy();
    const int group_i = region_i->getGroupIndex();

    if(region_i->getType() == RegionType::BLOCKAGE)
      continue;

    for(int j = i + 1; j < num_regions; j++)
    {
      auto& region_j = regions_[j];
      const int lx_j = region_j->getLx();
      const int ly_j = region_j->getLy();
      const int ux_j = region_j->getUx();
      const int uy_j = region_j->getUy();
      const int group_j = region_j->getGroupIndex();

      if(region_j->getType() == RegionType::BLOCKAGE)
        continue;

      if(group_i != group_j)
        continue;

      if(ux_i == lx_j or lx_i == ux_j or uy_i == ly_j or ly_i == uy_j)
      {
        region_i->addAdjacentRegion(region_j);
        region_j->addAdjacentRegion(region_i);
      }
    }
  }

  for(const auto region : regions_)
  {
    if(region->getType() == RegionType::BLOCKAGE)
      continue;

    if(region->isIsolated() == true)
    {
      const int lx = region->getLx();
      const int ly = region->getLy();
      const int ux = region->getUx();
      const int uy = region->getUy();
      //printf("Isolated (%6d, %6d) - (%6d, %6d)\n", lx, ly, ux, uy);
    }
  }
}

void
Legalizer::addNewRegion(std::shared_ptr<LGRegion> new_region)
{
  auto& regions_this_group = index2regions_[new_region->getGroupIndex()];
  regions_.push_back(new_region);
  regions_this_group.push_back(new_region);
}

void
Legalizer::handleRegionSegments()
{
  auto& pixels = grid_->getPixels();

  const int x_size = grid_->getXSize();
  const int y_size = grid_->getYSize();

  for(int i = 0; i < x_size; i++)
  {
    for(int j = 0; j < y_size; j++)
    {
      auto& pixel = pixels[i][j];
      if(pixel.h_overlap > 0 or pixel.v_overlap > 0)
      {
        if(pixel.h_overlap < site_width_ or pixel.v_overlap < row_height_)
          pixel.valid = false;
      }
    }
  }
}

} // namespace legalizer
