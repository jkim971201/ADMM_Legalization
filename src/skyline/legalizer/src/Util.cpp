#include <fstream>
#include <cstdio>
#include <thread>
#include <mutex>
#include <cassert>

#include "legalizer/Legalizer.h"

#include "Grid.h"
#include "objects/LGCell.h"
#include "objects/LGVioStat.h"
#include "objects/LGRegion.h"

#include "db/dbDatabase.h"
#include "db/dbTech.h"
#include "db/dbDesign.h"

namespace legalizer
{

inline std::vector<std::string>
tokenize(std::string_view line, std::string_view dels)
{
  std::string token;
  std::vector<std::string> tokens;

  for(auto itr = line.begin(); itr < line.end(); itr++)
  {
    bool is_del = dels.find(*itr) != std::string_view::npos;
    if(is_del || std::isspace(*itr))
    {
      if(!token.empty())
      {
        tokens.push_back(std::move(token));
        token.clear();
      }
    }
    else
      token.push_back(*itr);
  }

  if(!token.empty()) 
    tokens.push_back(std::move(token));

  // No need to move?
  // -> Return Value Optimization
  return tokens;
}

void
Legalizer::printDesignInfo() const
{
  const auto design = dbDatabase_->getDesign();
  const auto groups = design->getGroups();

  int num_cells          = static_cast<int>(cells_.size());
  int num_fixed_macros   = static_cast<int>(fixed_macro_rects_.size());
  int num_fixed_stdcells = static_cast<int>(fixed_stdcell_rects_.size());
  int num_groups         = static_cast<int>(groups.size());

  double sum_area = double(sum_cell_area_);
  double core_area = double(grid_->getArea());
  double fixed_area = double(sum_fixed_area_);
  double movable_area = double(sum_movable_area_);
  double white_space = core_area - fixed_area;
  //double density = movable_area / white_space * 100.0; // definition of density from OpenROAD
  double density = movable_area / white_space * 100.0; // definition of density from Legalization context
  double util = sum_area / core_area * 100.0;

  printf("------------------------------------------\n");
  printf("Design           : %13s\n", design->name().c_str());
  printf("# Cells          : %13d\n", num_cells);
  printf("# Fixed Macros   : %13d\n", num_fixed_macros);
  printf("# Fixed StdCells : %13d\n", num_fixed_stdcells);
  printf("# Groups         : %13d\n", num_groups);
  printf("# Total Rows     : %13d\n", num_total_rows_);
  printf("Site Width       : %13d\n", site_width_);
  printf("Row Height       : %13d\n", row_height_);
  printf("Utilization      : %13.1f %%\n", util);
  printf("Density          : %13.1f %%\n", density);
  if(!constraint_.empty())
  {
    printf("Max Util         : %13.1f %%\n", max_util_constraint_);
    printf("Max Disp         : %13.0f\n", max_disp_constraint_);
  } 

  //if(height_distribution_.size() > 1)
  {
    for(auto& [inst_height, num_cells_this_height] : height_distribution_)
    {
      const int height_multiplier = inst_height / row_height_;
      const double ratio = static_cast<double>(num_cells_this_height)
                         / static_cast<double>(num_cells) * 100.0;
      printf("%dH Cell          : %13d (%4.1f %%)\n", 
          height_multiplier, num_cells_this_height, ratio);
    }
  }

  int core_lx = grid_->getLx();
  int core_ly = grid_->getLy();
  int core_ux = grid_->getUx();
  int core_uy = grid_->getUy();
  printf("Core (%d, %d) - (%d, %d)\n", core_lx, core_ly, core_ux, core_uy);
  printf("------------------------------------------\n");
}

void
Legalizer::printStats(const LGStat& lg_stat, const LGVioStat& vio_stat, double rt) const
{
  const double dbu = static_cast<double>(dbDatabase_->getTech()->getDbu());
  const double sw_double = static_cast<double>(site_width_);
  const double rh_double = static_cast<double>(row_height_);
  const double delta_hpwl = lg_stat.hpwl_after / lg_stat.hpwl_before * 100.0 - 100.0;
  //const double delta_hpwl = (lg_stat.hpwl_after - lg_stat.hpwl_before) / lg_stat.hpwl_before * 100.0;
  printf("------------------------------------------\n");
  printf("Legalization Time : %12.1f s\n", rt);
  printf("\n");
  printf("Sum Displacement  : %12.1f um\n", lg_stat.total_disp);
  printf("Avg Displacement  : %12.1f um\n", lg_stat.avg_disp);
  printf("Max Displacement  : %12.1f um\n", lg_stat.max_disp);
  printf("\n");
  printf("HPWL Before       : %12.1f um\n", lg_stat.hpwl_before);
  printf("HPWL After        : %12.1f um\n", lg_stat.hpwl_after);
  printf("Delta HPWL        : %12.2f %%\n", delta_hpwl);
  printf("\n");
  printf("Avg Displacement  : %12.2f (rows) \n", lg_stat.avg_disp   * dbu / rh_double);
  printf("Avg Displacement  : %12.2f (sites)\n", lg_stat.avg_disp   * dbu / sw_double);
  printf("Max Displacement  : %12.2f (rows) \n", lg_stat.max_disp   * dbu / rh_double);
  printf("Max Displacement  : %12.2f (sites)\n", lg_stat.max_disp   * dbu / sw_double);
  printf("Sum Displacement  : %12.2f (sites)\n", lg_stat.total_disp * dbu / sw_double);

  printf("\n");
  printf("NumOverlap        : %12d\n", vio_stat.total_ovf);
  printf("NumRowOrientVio   : %12d\n", vio_stat.num_row_orient_vio);

  if(constraint_.empty() == false)
    printICCAD17Metric(delta_hpwl, vio_stat);
  printf("------------------------------------------\n");
}

void
Legalizer::printICCAD17Metric(double delta_hpwl_percentage, const LGVioStat& vio_stat) const
{
  std::map<int, double> h_to_sum_disp_sw;
  const double site_width_double = static_cast<double>(site_width_);
  const double row_height_double = static_cast<double>(row_height_);
  const double rh_sw_ratio = row_height_double / site_width_double;

  /* Compute S_am */
  int max_disp_dbu = 0;
  for(const auto& cell : cells_)
  {
    int cell_height = cell.getHeight();
    int disp_dbu_x = std::abs(cell.getLx() - cell.getInitLx());
    int disp_dbu_y = std::abs(cell.getLy() - cell.getInitLy());
    int disp_dbu = disp_dbu_x + disp_dbu_y;
    max_disp_dbu = std::max(max_disp_dbu, disp_dbu);

    double disp_in_sw = double(disp_dbu) / site_width_double;
    h_to_sum_disp_sw[cell_height] += disp_in_sw;
  }

  double total_disp_sw = 0.0;
  for(auto& [h, sum_disp_sw] : h_to_sum_disp_sw)
  {
    int num_k_cells = height_distribution_.at(h);
    total_disp_sw += sum_disp_sw / static_cast<double>(std::max(1, num_k_cells));
  }

  double large_k = static_cast<double>(height_distribution_.size());
  double S_am = total_disp_sw / large_k / rh_sw_ratio;
  double M_max = max_disp_dbu / row_height_double;

  /* Compute Total Score */
  int N_p = vio_stat.num_pin_short_vio + vio_stat.num_pin_access_vio;
  int N_e = vio_stat.num_edge_spacing_vio;
  double drv_penalty = static_cast<double>(N_p + N_e) / static_cast<double>(cells_.size());
  double S = (1.0 + delta_hpwl_percentage / 100.0 + drv_penalty) * (1.0 + M_max / 100.0) * S_am;

  printf("\n");
  printf("ICCAD 17 Contest Metric\n");
  printf("N_p   : %6d  \n", N_p);
  printf("N_e   : %6d  \n", N_e);
  printf("M_max : %6.3f\n", M_max);
  printf("S_am  : %6.3f\n", S_am);
  printf("S     : %6.3f\n", S);
}

void
Legalizer::parseICCAD17Constraint()
{
  std::string line = std::string();
  std::ifstream input_file(constraint_);

  while(std::getline(input_file, line))
  {
    const auto tokens = tokenize(line, "=");

    if(tokens.front() == "maximum_utilization")
    {
      std::string util_percent = tokens.back();
      std::string util_str = util_percent.substr(0, util_percent.size() - 1);
      // cut-off % character
      float util = std::stof(util_str);
      //printf("Util : %f\n", util);
      max_util_constraint_ = util;
    }
    if(tokens.front() == "maximum_movement")
    {
      std::string disp_rows_str = tokens.back();
      std::string disp_str = disp_rows_str.substr(0, disp_rows_str.size() - 4);
      float disp = std::stof(disp_str);
      //printf("Disp : %f\n", disp);
      max_disp_constraint_ = disp * static_cast<float>(row_height_);
    }
  }
}

void
Legalizer::parseSizeFile(std::unordered_map<std::string, std::pair<int, int>>& name_to_size)
{
  std::string line = std::string();
  std::ifstream input_file(size_file_);

  while(std::getline(input_file, line))
  {
    const auto tokens = tokenize(line, "=");
    std::string cell_name = tokens[0];
    int cell_width_in_site = std::stoi(tokens[1]);
    int cell_height_in_row = std::stoi(tokens[2]);
    name_to_size[cell_name] = {cell_width_in_site, cell_height_in_row};
  }
}

void
Legalizer::printRegionInfo() const
{
  std::ofstream output("log.txt");
  int x_size = grid_->getXSize();
  int y_size = grid_->getYSize();
  const auto& pixels = grid_->getPixels();

  for(int x = 0; x < x_size; x++)
    for(int y = 0; y < y_size; y++)
      if(pixels[x][y].region == nullptr)
        output << x << " " << y << std::endl;
  return;

  std::ofstream output2("log2.txt");
  for(auto region : regions_) 
  {
    if(region->getType() == RegionType::FENCE)
    {
      output2 << region->getLx() << " ";
      output2 << region->getLy() << " ";
      output2 << region->getWidth() << " ";
      output2 << region->getHeight() << std::endl;
    }
  }
}

void
Legalizer::printValidInfo() const
{
  const auto& pixels = grid_->getPixels();

  std::ofstream output("valid_map.txt");

  int x_size = grid_->getXSize();
  int y_size = grid_->getYSize();

  for(int x = 0; x < x_size; x++)
  {
    for(int y = 0; y < y_size; y++)
      output << x << " " << y << " " << int(pixels[x][y].valid) << std::endl;
  }
}

}
