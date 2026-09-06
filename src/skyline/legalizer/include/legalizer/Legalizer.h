#ifndef LEGALIZER_H 
#define LEGALIZER_H

#include <limits>
#include <vector>
#include <string>
#include <optional>
#include <memory>
#include <map>
#include <set>
#include <unordered_map>

namespace db {
  class dbDatabase;
  class dbTech;
  class dbLayer;
  class dbDesign;
  class dbTypes;
  class dbDie;
  class dbInst;
  class dbNet;
  class dbITerm;
  class dbBTerm;
  class dbMacro;
  class dbGroup;
  class dbRegion;
  class dbRow;
}

namespace legalizer 
{

constexpr int k_default_group_index = -1;
constexpr int k_blockage_group_index = -2;
constexpr int k_invalid_group_index = -3;

using namespace db;

class LGRect;
class LGCell;
class LGCellTechInfo;
class LGRegion;
class LGHyperParameters;

struct LGVioStat;

class Grid;
class CudaDatabase;

struct LGStat
{
  double max_disp = 0.0;
  double avg_disp = 0.0;
  double total_disp = 0.0;
  double hpwl_before;
  double hpwl_after;
};

struct ProjectionResult
{
  int cell_lx;
  int cell_ly;
  int disp;
};

class Legalizer 
{
  public:

    Legalizer();
    Legalizer(std::shared_ptr<db::dbDatabase> db);
    ~Legalizer();

    // APIs
    void run(); // Run legalization algorithm
    void setICCAD17Constraint(std::string_view file_path) { constraint_ = file_path; }
    void setSizeByFile(std::string_view file_path) { size_file_ = file_path; }
    void setTechPenalty(float val);
    void setStandardAdmm(bool val);
    void setQpRefine(bool val);
    void setCoeffToInitRho(float val);
    void setCoeffToInitDual(float val);
    void setMaxAdmmIter(int iter);

    void setFlagICCAD2017(bool val);

    void setSortingAlgorithm(std::string_view algo);
    void setSortingPolicy(std::string_view policy);

    void setAdaptiveSearch(bool val);
    void setMinIterToCallAdaptiveSearch(int iter);
    void setOvfRatioToCallAdaptiveSearch(float ratio);
    void setCoeffToNeglectDisplace(float val);

    void setRandomSorting(bool val);
    void setRandomSortingSeed(int seed);

  private:

    /* ------------------------ Methods ------------------------- */
    // Legalizer.cpp
    void initParameters();
    void adjustYHint();
    void makeDirectionVectors();
    void initializeCudaDatabase();
    void runAdmmLegalize();
    void runPostProcess();
    void detectIrregularShape();

    int64_t computeHpwlFromDB() const;
    LGStat computeStat() const;

    // DatabaseIO.cpp
    void importDB();
      void importCells();
      void importBlockageRegions();
      void importGroupRegions();
      void importNonGroupRegions();
      void importTechInfo();
        void importLayers(const std::vector<dbLayer*>& db_layers);
        void importEdgeSpcingRules(int dbu, 
          const std::vector<std::tuple<int, int, double>>& edge_spacing_rules);
        void importCellMacro(const std::vector<dbMacro*>& db_macros);
        void importPowerMetals(const std::vector<dbNet*>& special_nets);
      void assignGroupMembers();
      void handleRegionSegments();
      void findRegionAdjacency();
      void createGrid();
      void markRegion();
      void addNewRegion(std::shared_ptr<LGRegion> new_region);

    void exportDB() const;

    // Geometry.cpp
    bool isOutOfCore(const LGRect& rect) const;
    bool isOutOfRegion(const LGRect* region, const LGCell& cell, int lx, int ly) const;
    bool checkOverlap(const LGRect& rect1, const LGRect& rect2) const;
    bool checkRowViolation(const LGCell& cell) const;
    ProjectionResult findProjectionPoint(const LGCell& cell, const LGRect* region) const;

    // Preplace.cpp
    void initialLegalize();
    void preplaceCell(LGCell& cell);

    void roundToNearsetPoint(LGCell& cell);
    void moveOutsideUncorrectRegion(LGCell& cell, const LGRect* rect);

    bool checkHopeless(LGCell& cell);

    // Util.cpp
    void printDesignInfo() const;
    void printStats(const LGStat& lg_stat, const LGVioStat& vio_stat, double rt) const;
    void printICCAD17Metric(double delta_hpwl, const LGVioStat& vio_stat) const;
    void parseICCAD17Constraint();
    void parseSizeFile(std::unordered_map<std::string, std::pair<int, int>>& name_to_size);
    void printRegionInfo() const;
    void printValidInfo() const;

    /* ------------------------ Members ------------------------- */
    std::shared_ptr<db::dbDatabase> dbDatabase_;

    std::shared_ptr<LGHyperParameters> params_;

    std::vector<int> direction_default_;
    std::vector<int> direction_horizontal_;
    std::vector<int> direction_vertical_;
    std::vector<int> direction_wide_;

    std::vector<LGCell> cells_;               // Only movable standard cells
    std::vector<LGRect> fixed_macro_rects_;   // We only need shape & postion.
    std::vector<LGRect> fixed_stdcell_rects_; // We only need shape & postion.
    
    std::vector<std::shared_ptr<LGRegion>> regions_;
    std::vector<std::shared_ptr<LGCellTechInfo>> cell_tech_infos_;

    std::shared_ptr<Grid> grid_;
    std::unordered_map<const db::dbInst*, LGCell*> db_inst_to_lg_cell_;
    std::unordered_map<const db::dbGroup*, int> db_group_to_index_;
    std::unordered_map<const db::dbLayer*, int> db_layer_to_index_;
    std::unordered_map<const db::dbMacro*, std::shared_ptr<LGCellTechInfo>> db_macro_to_tech_info_;
    std::unordered_map<int, std::vector<std::shared_ptr<LGRegion>>> index2regions_;

    std::shared_ptr<CudaDatabase> cuda_database_;

    // Design Info
    int core_lx_;
    int core_ly_;
    int core_ux_;
    int core_uy_; 
    int site_width_;           // Assume Uniform SiteWidth
    int row_height_;           // Assume Uniform RowHeight
    int num_total_rows_;       // Num of Total Rows 
    int64_t original_hpwl_;    // Initial HPWL
    int64_t sum_cell_area_;    // Sum Area of Total Instance
    int64_t sum_movable_area_; // Sum Area of Movable Cells
    int64_t sum_fixed_area_;   // Sum Area of Fixed Macros

    std::map<int, int> height_distribution_; // (key, val) : (inst height, #insts)

    std::vector<int> edge_spacing_rules_; 
    // 0 : Empty
    // 1 : Empty
    // 1 + 1 : Spacing (in site) between type1 and type1
    // 1 + 2 : Spacing (in site) between type1 and type2
    // 2 + 2 : Spacing (in site) between type2 and type2
    // 5:  Empty
    // <- I know this is stupid, but I cannot think of better idea for array-based representation.
    // (since CUDA kernel only accepts array-style input)

    // Size Info (only used by modified ISPD 2015 Benchmarks)
    std::string size_file_; 

    // Constraint Info (can be empty unless using ICCAD17 Benchmarks)
    std::string constraint_; // .constraint file of ICCAD17 Contest Benchmark
    double max_disp_constraint_;
    double max_util_constraint_;
    /* ---------------------------------------------------------- */
};

} // namespace legalizer 

#endif
