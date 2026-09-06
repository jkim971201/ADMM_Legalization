#include <cmath>
#include <cassert>
#include <map>
#include <random>    // For mt19937
#include <algorithm> // For sort
#include <cstdio>
#include <iostream>
#include <fstream>

#include "db/dbDatabase.h"
#include "db/dbTech.h"
#include "db/dbDesign.h"
#include "db/dbTypes.h"
#include "db/dbDie.h"
#include "db/dbInst.h"
#include "db/dbNet.h"
#include "db/dbITerm.h"
#include "db/dbBTerm.h"

#include "object/GPObject.h"
#include "SkyPlaceDB.h"
#include "Util.h"

namespace skyplace
{

inline float getOverlapArea(const GPBin* bin, const GPCell* cell)
{
  float rectLx = std::max(bin->lx(), cell->lx());
  float rectLy = std::max(bin->ly(), cell->ly());
  float rectUx = std::min(bin->ux(), cell->ux());
  float rectUy = std::min(bin->uy(), cell->uy());

  if(rectLx >= rectUx || rectLy >= rectUy)
    return 0;
  else
    return (rectUx - rectLx) * (rectUy - rectLy);
}

// getOverlapArea should use int64_t ideally,
// but runtime will be doubled (referred to OpenROAD comment)
inline float getOverlapAreaWithDensitySize(const GPBin* bin, const GPCell* cell)
{
  float rectLx = std::max(bin->lx(), cell->dLx());
  float rectLy = std::max(bin->ly(), cell->dLy());
  float rectUx = std::min(bin->ux(), cell->dUx());
  float rectUy = std::min(bin->uy(), cell->dUy());

  if(rectLx >= rectUx || rectLy >= rectUy)
    return 0;
  else
    return ((rectUx - rectLx) * (rectUy - rectLy));
}

SkyPlaceDB::SkyPlaceDB()
  : designName_          (""),
    targetDensity_       (1.0),
    dbu_                 (0),
    numStdCells_         (0), 
    numMacro_            (0),
    numFixed_            (0), 
    numMovable_          (0), 
    numMovableMacro_     (0),
    numFiller_           (0),
    numBinX_             (0), 
    numBinY_             (0), 
    numIO_               (0),
    numCluster_          (0),
    sumFixedArea_        (0),
    sumMovableArea_      (0), 
    sumMovableStdArea_   (0), 
    sumMovableMacroArea_ (0), 
    sumScaledMovableArea_(0), 
    sumTotalInstArea_    (0),
    dieArea_             (0),
    density_             (0),
    util_                (0),
    avgCellArea_         (0),
    avgStdCellArea_      (0),
    hpwl_                (0), 
    binX_                (0), 
    binY_                (0),
    die_ptr_             (nullptr),
    fillerWidth_         (0), 
    fillerHeight_        (0),
    numInitStep_         (5),
    is_bookshelf_        (false) 
{
  reset();
}

void
SkyPlaceDB::reset()
{
  dbu_ = 0;
  
  designName_ = std::string();

  die_ptr_ = nullptr;

  cellPtrs_.clear();
  cellInsts_.clear();

  movableMacroPtrs_.clear();
  fixedPtrs_.clear();
  movablePtrs_.clear();

  netPtrs_.clear();
  netInsts_.clear();

  pinPtrs_.clear();
  pinInsts_.clear();

  binPtrs_.clear();
  binInsts_.clear();

  dbInst2Cell_.clear();
}

void
SkyPlaceDB::exportDB(std::shared_ptr<dbDatabase> _dbDatabase)
{
  // Due to the bookshelf conversion
  const float scaling_factor = static_cast<float>(bookshelf_dbu_); 
  for(auto cell : movablePtrs_)
  {
    if(!cell->isFiller())
    {
      dbInst* inst_ptr = cell->dbInstPtr();
      float new_lx_float = std::round(cell->lx() * scaling_factor);
      float new_ly_float = std::round(cell->ly() * scaling_factor);

      int new_lx_int = static_cast<int>(new_lx_float);
      int new_ly_int = static_cast<int>(new_ly_float);
      inst_ptr->setLocation(new_lx_int, new_ly_int);
    }
  }
}

void
SkyPlaceDB::importDB(std::shared_ptr<dbDatabase> _dbDatabase)
{
  const auto tech   = _dbDatabase->getTech();
  const auto design = _dbDatabase->getDesign();

  dbu_           = tech->getDbu();
  designName_    = design->name();
  bookshelf_dbu_ = _dbDatabase->getBookshelfDbu();
  is_bookshelf_  = _dbDatabase->isBookShelf();

  const std::vector<dbNet*>&   db_nets   = design->getNets();
  const std::vector<dbInst*>&  db_insts  = design->getInsts();
  const std::vector<dbITerm*>& db_iterms = design->getITerms();
  const std::vector<dbBTerm*>& db_bterms = design->getBTerms();

  const float scaling_factor = static_cast<float>(bookshelf_dbu_); 

  // Step1. Initialize Die
  die_ptr_ = std::make_shared<GPDie>(
    design->coreLx() / scaling_factor,
    design->coreLy() / scaling_factor,
    design->coreUx() / scaling_factor,
    design->coreUy() / scaling_factor);

  // Step2. Initialize GPCell
  int numInst = db_insts.size();

  int fixedIdx   = 0;
  int movableIdx = 0;

  cellInsts_.reserve(numInst);
  cellPtrs_.reserve(numInst);

  // values will be divided by this value
  for(auto& db_inst_ptr : db_insts)
  {
    cellInsts_.emplace_back(db_inst_ptr, scaling_factor);
    GPCell* new_cell_ptr = &cellInsts_.back();

    dbInst2Cell_[db_inst_ptr] = new_cell_ptr;
    cellPtrs_.push_back(new_cell_ptr);

    if(new_cell_ptr->isFixed() == true) 
    {
      new_cell_ptr->setID(fixedIdx++);
      fixedPtrs_.push_back(new_cell_ptr);
      if(!isOutsideDie(new_cell_ptr))
        sumFixedArea_ += new_cell_ptr->area();
      if(new_cell_ptr->isMacro() == true)
        numMacro_++;
    }
    else                    
    {
      new_cell_ptr->setID(movableIdx++);
      movablePtrs_.push_back(new_cell_ptr);
      auto cellArea = new_cell_ptr->area();
      sumMovableArea_ += cellArea;

      if(new_cell_ptr->isMacro() == true)
      {
        sumMovableMacroArea_  += cellArea;
        sumScaledMovableArea_ += cellArea * targetDensity_;
        numMacro_++;
        numMovableMacro_++;
      }
      else 
      {
        sumMovableStdArea_    += cellArea;
        sumScaledMovableArea_ += cellArea;
        numStdCells_++;
      }
    }
  }

  // Step3. Initialize GPNet
  int numNets = db_nets.size();
  netInsts_.reserve(numNets);
  netPtrs_.reserve(numNets);

  for(auto& db_net_ptr : db_nets)
  {
    if(db_net_ptr->isSpecial() == true)
      continue;

    netInsts_.emplace_back(db_net_ptr, netInsts_.size());
    GPNet* new_net_ptr = &netInsts_.back();
    netPtrs_.push_back(new_net_ptr);
  }

  // Step4. Initialize GPPin
  const int num_bterms = static_cast<int>(db_bterms.size());
  const int num_iterms = static_cast<int>(db_iterms.size());
  const int num_pins   = num_bterms + num_iterms;
  pinInsts_.reserve(num_pins);
  pinPtrs_.reserve(num_pins);
  io_pins_.reserve(num_bterms);

  auto cell_table_end = dbInst2Cell_.end();
  for(auto& net_ptr : netPtrs_)
  {
    dbNet* db_net_ptr = net_ptr->dbNetPtr();
    if(db_net_ptr->isSpecial() == true)
      continue;

    auto& bterms_this_net = db_net_ptr->getBTerms();
    for(auto& bterm : bterms_this_net)
    {
      pinInsts_.emplace_back(bterm, pinInsts_.size(), scaling_factor);
  
      GPPin* new_pin_ptr = &pinInsts_.back();
      pinPtrs_.push_back(new_pin_ptr);
      io_pins_.push_back(new_pin_ptr);

      net_ptr->addNewPin(new_pin_ptr);
      new_pin_ptr->setNet(net_ptr);
    }
  
    auto& iterms_this_net = db_net_ptr->getITerms();
    for(auto& iterm : iterms_this_net)
    {
      pinInsts_.emplace_back(iterm, pinInsts_.size(), scaling_factor);
  
      GPPin* new_pin_ptr = &pinInsts_.back();
      pinPtrs_.push_back(new_pin_ptr);

      net_ptr->addNewPin(new_pin_ptr);
      new_pin_ptr->setNet(net_ptr);

      auto cell_find = dbInst2Cell_.find(iterm->getInst());
      GPCell* cell_ptr = cell_find->second;

      cell_ptr->addNewPin(new_pin_ptr);
      new_pin_ptr->setCell(cell_ptr);
    }
  }

  // Step5. Initialize GPPin Location and NetBBox
  updateHpwl();

  // Finish -> assign number variables
  numMovable_ = static_cast<int>(movablePtrs_.size());
  numFixed_   = static_cast<int>(fixedPtrs_.size());
  numIO_      = static_cast<int>(db_bterms.size());

  sumTotalInstArea_ = sumFixedArea_ + sumMovableArea_;
  dieArea_ = die_ptr_->area();
  density_ = sumTotalInstArea_ / dieArea_;
  util_    = sumMovableArea_   / (dieArea_ - sumFixedArea_);

  avgCellArea_    = (sumMovableArea_ + sumFixedArea_) / (numMovable_  + numFixed_);
  avgStdCellArea_ = (sumMovableArea_) / static_cast<float>(numMovable_);

  // If target density is smaller than util, fix it.
  // NOTE : This must be done before createBins()
  validateTargetDensity();
}

void
SkyPlaceDB::validateTargetDensity()
{
  float util_percentage                   = util_ * 100.0;
  float target_density_percentage         = targetDensity_ * 100.0;
  float target_density_default            = 1.0;
  float target_density_default_percentage = target_density_default * 100.0;

  // In the bookshelf flow, there are some cases that
  // target density is lower than util (ex. newblue4)
  if(is_bookshelf_ == false) 
  {
    if(target_density_percentage <= util_percentage)
    {
      printf("Target density (%4.2f%%) is less than Util (%4.2f%%).\n", 
          target_density_percentage, util_percentage);
      exit(0);
    }
  }
}

void
SkyPlaceDB::init(std::shared_ptr<dbDatabase> db)
{
  auto db_init_start = std::chrono::high_resolution_clock::now();

  reset();

  // Step#1: Import dbDatabase to SkyPlaceDB
  importDB(db);
  const double db_import_time = evalTime(db_init_start);
  printf("Import     Database           (takes %5.2f s)\n", db_import_time);

  // Step#2: GPBin Initialization
  auto create_bin_start = std::chrono::high_resolution_clock::now();
  createBins();
  const double bin_init_time = evalTime(create_bin_start);
  printf("Initialize BinGrid            (takes %5.2f s)\n", bin_init_time);
  
  // Step#3: Filler Insertion
  auto create_filler_start = std::chrono::high_resolution_clock::now();
  createFillers();
  const double filler_time = evalTime(create_filler_start);
  printf("Initialize FillerCell         (takes %5.2f s)\n", filler_time);

  // Step#4: Update Density Size
  auto density_size_start = std::chrono::high_resolution_clock::now();
  updateDensitySize();
  const double density_size_time = evalTime(density_size_start);
  printf("Initialize DensitySize        (takes %5.2f s)\n", density_size_time);

  // Step#5: Update Fixed Overlap Area
  auto fixed_overlap_start = std::chrono::high_resolution_clock::now();
  updateFixedOverlapArea();
  const double fixed_overlap_time = evalTime(fixed_overlap_start);
  printf("Initialize FixedOverlap       (takes %5.2f s)\n", fixed_overlap_time);
}

void
SkyPlaceDB::updateHpwl()
{
  // Initialize GPPin Location
  for(auto& cell_ptr : cellPtrs_)
    for(auto& pin_ptr : cell_ptr->pins())
      pin_ptr->updatePinLocation(cell_ptr);

  hpwl_ = 0;

  // Initilize GPNet BBox
  for(auto& net_ptr : netPtrs_)
  {
    net_ptr->updateBBox();
    hpwl_ += net_ptr->hpwl();
  }
}

bool
SkyPlaceDB::isOutsideDie(const GPCell* cell)
{
  bool isOutside = false;

  if(cell->ux() <= die_ptr_->lx()) isOutside = true;
  if(cell->lx() >= die_ptr_->ux()) isOutside = true;
  if(cell->uy() <= die_ptr_->ly()) isOutside = true;
  if(cell->ly() >= die_ptr_->uy()) isOutside = true;

  return isOutside;
}

void
SkyPlaceDB::updatePinBound()
{
  for(auto& net : netPtrs_)
  {
    GPPin* minPinX = nullptr;
    GPPin* minPinY = nullptr;
    GPPin* maxPinX = nullptr;
    GPPin* maxPinY = nullptr;

    float minX = std::numeric_limits<float>::max();
    float minY = std::numeric_limits<float>::max();
    float maxX = 0;
    float maxY = 0;

    for(auto& pin : net->pins())
    {
      float cx = pin->cx();
      float cy = pin->cy();

      if(cx <= minX)
      {
        if(minPinX != nullptr)
          minPinX->unsetMinPinX();
        minX = cx;
        pin->setMinPinX();
        minPinX = pin;
      }

      if(cy <= minY)
      {
        if(minPinY != nullptr)
          minPinY->unsetMinPinY();
        minY = cy;
        pin->setMinPinY();
        minPinY = pin;
      }

      if(cx > maxX)
      {
        if(maxPinX != nullptr)
          maxPinX->unsetMaxPinX();
        maxX = cx;
        pin->setMaxPinX();
        maxPinX = pin;
      }

      if(cy > maxY)
      {
        if(maxPinY != nullptr)
          maxPinY->unsetMaxPinY();
        maxY = cy;
        pin->setMaxPinY();
        maxPinY = pin;
      }
    }
  }
}

void
SkyPlaceDB::createBins()
{
  // Do not use density_, rather use TargetDensity 
  idealBinAreaForAvgCellArea_    = avgCellArea_    / targetDensity_;
  idealBinAreaForAvgStdCellArea_ = avgStdCellArea_ / targetDensity_;

  idealBinCountForAvgCellArea_    = dieArea_ / idealBinAreaForAvgCellArea_;
  idealBinCountForAvgStdCellArea_ = dieArea_ / idealBinAreaForAvgStdCellArea_;

  idealNumBinForAvgCellArea_      = std::sqrt(idealBinCountForAvgCellArea_);
  idealNumBinForAvgStdCellArea_   = std::sqrt(idealBinCountForAvgStdCellArea_);

  //idealNumBin_ = std::ceil(idealNumBinForAvgCellArea_);
  idealNumBin_ = std::ceil(idealNumBinForAvgStdCellArea_);

  float dieX = die_ptr_->dx();
  float dieY = die_ptr_->dy();

  float aspect_ratio = std::max(dieX / dieY, dieY / dieX);
  int ratio = 1;

  if(aspect_ratio >= 4.0)
    ratio = 4;
  else if(aspect_ratio >= 2.0)
    ratio = 2;
  else 
    ratio = 1;

  int numBin = 4; 
  while(numBin * numBin * ratio < idealNumBin_ * idealNumBin_)
  {
    numBin *= 2;
    if(numBin >= 1024)
      break;
  }

  int error_num_bin1 
    = std::abs(idealNumBin_ * idealNumBin_ - numBin * numBin * ratio);
  int error_num_bin2 
    = std::abs(idealNumBin_ * idealNumBin_ - numBin * numBin * ratio / 4);
  if(error_num_bin2 < error_num_bin1)
    numBin = numBin / 2;

  if(dieX > dieY)
  {
    numBinX_ = numBin * ratio;
    numBinY_ = numBin;
  }
  else
  {
    numBinX_ = numBin;
    numBinY_ = numBin * ratio;
  }

  binInsts_.resize(numBinX_ * numBinY_);
  binPtrs_.reserve(numBinX_ * numBinY_);

  float binSizeX = die_ptr_->dx() / static_cast<float>(numBinX_);
  float binSizeY = die_ptr_->dy() / static_cast<float>(numBinY_);

  float dieUx = die_ptr_->ux();
  float dieLx = die_ptr_->lx();
  float dieUy = die_ptr_->uy();
  float dieLy = die_ptr_->ly();

  float lx = dieLx;
  float ly = dieLy;
  int numCreated = 0;

  for(auto& bin : binInsts_)
  {  
    float binWidth  = binSizeX;
    float binHeight = binSizeY;

    int row = numCreated / numBinX_;
    int col = numCreated % numBinX_;

    if(col == (numBinX_ - 1) || (lx + binSizeX > dieUx))
      binWidth  = dieUx - lx;

    if(row == (numBinY_ - 1) || (ly + binSizeY > dieUy))
      binHeight = dieUy - ly;

    bin = GPBin(row, col, 
                lx, ly, 
                lx + binWidth, ly + binHeight, 
                targetDensity_);

    lx += binWidth;
    
    if(lx >= dieUx)
    {
      ly += binHeight;
      lx = dieLx;
    }

    binPtrs_.push_back(&bin);
    numCreated++;
  }

  binX_  = binSizeX;
  binY_  = binSizeY;
}

void
SkyPlaceDB::createFillers()
{
  float dxSum = 0;
  float dySum = 0;

  std::vector<float> dxList;
  std::vector<float> dyList;

  dxList.reserve(numMovable_);
  dyList.reserve(numMovable_);

  for(auto& cell : movablePtrs_)
  {
    dxList.push_back(cell->dx());
    dyList.push_back(cell->dy());
  }

  std::sort(dxList.begin(), dxList.end());

  int minIdx = static_cast<int>(static_cast<float>(dxList.size()) * 0.10);
  int maxIdx = static_cast<int>(static_cast<float>(dxList.size()) * 0.90);

  // if numMovable is too small
  if(minIdx == maxIdx) 
  {
    minIdx = 0;
    maxIdx = numMovable_;
  }

  for(int i = minIdx; i < maxIdx; i++)
  {
    dxSum += dxList[i];
    dySum += dyList[i];
  }

  fillerWidth_   = dxSum / static_cast<float>(maxIdx - minIdx);
  fillerHeight_  = dySum / static_cast<float>(maxIdx - minIdx);

  avgFillerArea_  = fillerWidth_ * fillerHeight_;
  whiteSpaceArea_ = (die_ptr_->area() - sumFixedArea_);
  fillerArea_     = whiteSpaceArea_ * targetDensity_ 
                  - sumScaledMovableArea_;

  if(fillerArea_ < 0)
  {
    printf("Error - FillerArea is smaller than 0!\n");
    printf("Use higher placement density...\n");
    exit(0);
  }

  numFiller_ = static_cast<int>(fillerArea_ / avgFillerArea_);

  std::mt19937 randVal(0);

  int dieLx = static_cast<int>(die_ptr_->lx());
  int dieLy = static_cast<int>(die_ptr_->ly());
  int dieWidth  = static_cast<int>(die_ptr_->dx());
  int dieHeight = static_cast<int>(die_ptr_->dy());

  // Rest of this function is incredibly stupid...
  // TODO: Fix these stupid parts
  for(int i = 0; i < numFiller_; i++)
  {
    auto randX = randVal();
    auto randY = randVal();
    
    // Random distribution over the entire layout
    GPCell filler(randX % dieWidth  + dieLx, 
                randY % dieHeight + dieLy,
                fillerWidth_, fillerHeight_);

    cellInsts_.push_back(filler);
  }

  cellPtrs_.clear();
  cellPtrs_.reserve(cellInsts_.size());

  numMovable_ += numFiller_;
  
  movableMacroPtrs_.clear();
  movableMacroPtrs_.reserve(numMovableMacro_);

  movablePtrs_.clear();
  movablePtrs_.reserve(numMovable_);

  fixedPtrs_.clear();
  fixedPtrs_.reserve(numFixed_);

  int fixedID   = 0;
  int movableID = 0;

  for(auto& cell : cellInsts_)
  {
    cellPtrs_.push_back(&cell);
    if(!cell.isFixed())
    {
      cell.setID(movableID++);
      movablePtrs_.push_back(&cell);

      if( cell.isMacro() )
        movableMacroPtrs_.push_back(&cell);
    }
    else
    {
      cell.setID(fixedID++);
      fixedPtrs_.push_back(&cell);
    }
    for(auto & pin : cell.pins())
      pin->setCell(&cell);
  }
}

OverlapBins
SkyPlaceDB::findBin(const GPCell* cell)
{
  float lx = cell->lx();
  float ux = cell->ux();
  
  int minX = std::floor((lx - die_ptr_->lx()) / binX_);
  int maxX = std::ceil((ux - die_ptr_->lx()) / binX_);

  minX = std::max(minX, 0);
  maxX = std::min(numBinX_, maxX);

  std::pair<int, int> minMaxX = std::make_pair(minX, maxX);

  float ly = cell->ly();
  float uy = cell->uy();
  
  int minY = std::floor((ly - die_ptr_->ly()) / binY_);
  int maxY = std::ceil((uy - die_ptr_->ly()) / binY_);

  minY = std::max(minY, 0);
  maxY = std::min(numBinY_, maxY);

  std::pair<int, int> minMaxY = std::make_pair(minY, maxY);

  return std::make_pair(minMaxX, minMaxY);
}

OverlapBins
SkyPlaceDB::findBinWithDensitySize(const GPCell* cell)
{
  float lx = cell->dLx();
  float ux = cell->dUx();
  
  int minX = std::floor((lx - die_ptr_->lx()) / binX_);
  int maxX = std::ceil((ux - die_ptr_->lx()) / binX_);

  minX = std::max(minX, 0);
  maxX = std::min(numBinX_, maxX);

  std::pair<int, int> minMaxX = std::make_pair(minX, maxX);

  float ly = cell->dLy();
  float uy = cell->dUy();

  int minY = std::floor((ly - die_ptr_->ly()) / binY_);
  int maxY = std::ceil((uy - die_ptr_->ly()) / binY_);

  minY = std::max(minY, 0);
  maxY = std::min(numBinY_, maxY);

  std::pair<int, int> minMaxY = std::make_pair(minY, maxY);

  return std::make_pair(minMaxX, minMaxY);
}

bool
SkyPlaceDB::isOverlap(const GPCell* cell, const GPBin* bin)
{
  bool checkX = false;
  bool checkY = false;

  if((cell->lx() >= bin->lx()) && (cell->lx() <= bin->ux()))
    checkX = true;
  if((cell->ux() >= bin->lx()) && (cell->ux() <= bin->ux()))
    checkX = true;
  if((cell->ly() >= bin->ly()) && (cell->ly() <= bin->uy()))
    checkY = true;
  if((cell->uy() >= bin->ly()) && (cell->uy() <= bin->uy()))
    checkY = true;

  return checkX && checkY;
}

void
SkyPlaceDB::moveCellInsideLayout(GPCell* cell)
{
  cell->setCenterLocation(getXCoordiInsideLayout(cell), 
                          getYCoordiInsideLayout(cell));
}

float
SkyPlaceDB::getXCoordiInsideLayout(const GPCell* cell)
{
  float newCx = cell->cx();

  if(cell->lx() < die_ptr_->lx())
    newCx = die_ptr_->lx() + cell->dx()/2;
  if(cell->ux() > die_ptr_->ux())
    newCx = die_ptr_->ux() - cell->dx()/2;

  return newCx;
}

float
SkyPlaceDB::getYCoordiInsideLayout(const GPCell* cell)
{
  float newCy = cell->cy();

  if(cell->ly() < die_ptr_->ly())
    newCy = die_ptr_->ly() + cell->dy()/2;
  if(cell->uy() > die_ptr_->uy())
    newCy = die_ptr_->uy() - cell->dy()/2;

  return newCy;
}

float
SkyPlaceDB::getXCoordiInsideLayout(const GPCell* cell, float cx)
{
  float newCx = cx;

  if(cx - cell->dx()/2 < die_ptr_->lx())
    newCx = die_ptr_->lx() + cell->dx()/2;
  if(cx + cell->dx()/2 > die_ptr_->ux())
    newCx = die_ptr_->ux() - cell->dx()/2;
  return newCx;
}

float
SkyPlaceDB::getYCoordiInsideLayout(const GPCell* cell, float cy)
{
  float newCy = cy;

  if(cy - cell->dy()/2 < die_ptr_->ly())
    newCy = die_ptr_->ly() + cell->dy()/2;
  if(cy + cell->dy()/2 > die_ptr_->uy())
    newCy = die_ptr_->uy() - cell->dy()/2;
  return newCy;
}

float
SkyPlaceDB::getXDensityCoordiInsideLayout(const GPCell* cell)
{
  float newCx = cell->cx();

  if(cell->dLx() < die_ptr_->lx())
    newCx = die_ptr_->lx() + cell->dDx()/2;
  if(cell->dUx() > die_ptr_->ux())
    newCx = die_ptr_->ux() - cell->dDx()/2;

  return newCx;
}

float
SkyPlaceDB::getYDensityCoordiInsideLayout(const GPCell* cell)
{
  float newCy = cell->cy();

  if(cell->dLy() < die_ptr_->ly())
    newCy = die_ptr_->ly() + cell->dDy()/2;
  if(cell->dUy() > die_ptr_->uy())
    newCy = die_ptr_->uy() - cell->dDy()/2;

  return newCy;
}

void
SkyPlaceDB::updateDensitySize()
{
  float scaleX = 0, scaleY = 0;
  float densityW = 0, densityH = 0;

  for(auto& cell: cellPtrs_)
  {
    if(cell->dx() > DENSITY_SCALE * binX_)
    {
      scaleX   = 1.0;
      densityW = cell->dx();
    }
    else
    {
      scaleX   = cell->dx() / (DENSITY_SCALE * binX_);
      densityW = DENSITY_SCALE * binX_;
    }

    if(cell->dy() > DENSITY_SCALE * binY_)
    {
      scaleY   = 1.0;
      densityH = cell->dy();
    }
    else
    {
      scaleY   = cell->dy() / (DENSITY_SCALE * binY_);
      densityH = DENSITY_SCALE * binY_;
    }
    cell->setDensitySize(densityW, densityH, scaleX * scaleY);
  }
}

void
SkyPlaceDB::updateMacroDensityWeight(float macroWeight)
{
  if(numMovableMacro_ == 0)
    return;

  float avgMacroArea = sumMovableMacroArea_ 
                     / static_cast<float>(numMovableMacro_);

  float maxMacroArea = 0.0;
  float minMacroArea = std::numeric_limits<float>::max();

  float sumMacroAreaVarianceSquare = 0.0;

  for(auto& macro : movableMacroPtrs_)
  {
    float area = macro->area();
    sumMacroAreaVarianceSquare += (area - avgMacroArea) 
                                * (area - avgMacroArea);
  
    if(area > maxMacroArea)
      maxMacroArea = area;
    if(area < minMacroArea)
      minMacroArea = area;
  }

  float std_dev = std::sqrt( sumMacroAreaVarianceSquare / 
         static_cast<float>( numMovableMacro_ ) );

  float maxRatio = maxMacroArea / avgMacroArea;
  float minRatio = minMacroArea / avgMacroArea;

  float ratio_std_macro = sumMovableMacroArea_ / sumMovableStdArea_;

  std::cout << "Avg MacroArea   : " << avgMacroArea << std::endl;
  std::cout << "Std_Dev         : " << std_dev << std::endl;
  std::cout << "Max Ratio       : " << maxRatio << std::endl;
  std::cout << "Min Ratio       : " << minRatio << std::endl;
  std::cout << "A_Macro / A_Std : " << ratio_std_macro << std::endl;
  std::cout << "Density Weight  : " << macroWeight << std::endl;
  std::cout << "Max MacroArea   : " << maxMacroArea << std::endl;
  std::cout << "Min MacroArea   : " << minMacroArea << std::endl;

  // if(targetDensity_ < 1.0 && (maxRatio > 1.2 || minRatio < 0.2) )
  for(auto& macro : movableMacroPtrs_)
  {
    //if(macro->dx() > die_ptr_->dx() * 0.2 || macro->dy() > die_ptr_->dy() * 0.2)
    {
      //printf("large Macro!\n");
      float weight = macro->densityScale();
      macro->setDensityScale(weight * macroWeight); // 1.09??
    }
  }
}

void
SkyPlaceDB::updateFixedOverlapArea()
{
  for(auto& bin : binPtrs_)
    bin->setFixedArea(0);

  for(auto& cell : fixedPtrs_) 
  {
    if(isOutsideDie(cell) == true)
      continue;

    OverlapBins ovBins = findBin(cell);
    int minX  = ovBins.first.first;
    int maxX  = ovBins.first.second;
    int minY  = ovBins.second.first;
    int maxY  = ovBins.second.second;

    for(int i = minX; i < maxX; i++)
    {
      for(int j = minY; j < maxY; j++)
      {
        GPBin* bin = binPtrs_[j * numBinX_ + i];
        const float overlapArea = getOverlapArea(bin, cell);
        bin->addFixedArea(overlapArea * bin->targetDensity());
        // FixedArea should be scaled-down with target density
        // (according to the OpenROAD RePlAce comment)
      }
    }
  }
}

void
SkyPlaceDB::printInfo() const 
{
  using namespace std;

  float initHpwl = hpwl_ / static_cast<float>(dbu_);
  // Zero Check is done in init()

  cout << endl;
  cout << "*** Summary of SkyPlaceDB ***" << endl;
  cout << "---------------------------------------------" << endl;
  cout << " DESIGN NAME        : " << designName_         << endl;
  cout << " NUM CELL (TOTAL)   : " << numMovable_ + numFixed_ << endl;
  cout << " NUM CELL (MOVABLE) : " << numMovable_         << endl;
  cout << " NUM CELL (FIXED)   : " << numFixed_           << endl;
  cout << " NUM CELL (FILLER)  : " << numFiller_          << endl;
  cout << " NUM MOVABLE MACRO  : " << numMovableMacro_    << endl;
  cout << " NUM NET            : " << numNet()            << endl;
  cout << " NUM PIN            : " << numPin()            << endl;
  cout << " NUM IO             : " << numIO()             << endl;
  cout << " NUM BIN   (IDEAL)  : " << idealNumBin_        << endl;
  cout << " NUM BIN   (TOTAL)  : " << numBinX_ * numBinY_ << endl;
  cout << " NUM BIN_X (USED)   : " << numBinX_            << endl;
  cout << " NUM BIN_Y (USED)   : " << numBinY_            << endl;
  cout << " BIN WIDTH          : " << binX_               << endl;
  cout << " BIN HEIGHT         : " << binY_               << endl;
  cout << " FILLER WIDTH       : " << fillerWidth_        << endl;
  cout << " FILLER HEIGHT      : " << fillerHeight_       << endl;
  cout << " FILLER AREA        : " << avgFillerArea_      << endl;
  cout << " AREA (TOTAL)       : " << sumTotalInstArea_   << endl;
  cout << " AREA (MOVABLE)     : " << sumMovableArea_     << endl;
  cout << " AREA (STD)         : " << sumMovableStdArea_  << endl;
  cout << " AREA (MACRO)       : " << sumMovableMacroArea_<< endl;
  cout << " AREA (FIXED)       : " << sumFixedArea_       << endl;
  cout << " AREA (FILLER)      : " << fillerArea_         << endl;
  cout << " AREA (WHITESPACE)  : " << whiteSpaceArea_     << endl;
  cout << " AREA (CORE)        : " << dieArea_            << endl;
  cout << " TARGET DENSITY     : " << targetDensity_ * 100.0 << "%\n";
  cout << " DENSITY            : " << density_       * 100.0 << "%\n";
  cout << " UTIL               : " << util_          * 100.0 << "%\n";
  cout << " INITIAL HPWL       : " << initHpwl << endl;
  cout << "---------------------------------------------" << endl;
  cout << endl;
}

void
SkyPlaceDB::printFixedOverlapArea()
{
  std::ofstream output;
  output.open("fixedOverlapArea.txt");

  for(auto& bin : binPtrs_)
    output << bin->fixedArea() << std::endl;

  output.close();
}

void
SkyPlaceDB::printBinDensity()
{
  std::ofstream output;
  output.open("BinDensity.txt");

  for(auto& bin : binPtrs_)
  {
    output << bin->col() << "," << bin->row() << ",";
    output << bin->density() << std::endl;
  }

  output.close();
}

void
SkyPlaceDB::printBinPotential()
{
  std::ofstream output;
  output.open("Potential.txt");

  for(auto& bin : binPtrs_)
  {
    output << bin->col() << "," << bin->row() << ",";
    output << bin->potential() << std::endl;
  }

  output.close();
}

void
SkyPlaceDB::printBinElectroFieldXY()
{
  std::ofstream output;
  output.open("ElectroField.txt");

  for(auto& bin : binPtrs_)
  {
    output << bin->col() << "," << bin->row() << ",";
    output << bin->electroForceX() << ",";
    output << bin->electroForceY() <<  std::endl;
  }

  output.close();
}

void
SkyPlaceDB::debugCell() const
{
  printf("=== Debug GPCell ===\n");
  printf("NumCells : %ld\n", cellPtrs_.size());
  for(const auto& cell : cellPtrs_)
    printCellInfo(cell);
}

void
SkyPlaceDB::debugNet() const
{
  printf("=== Debug GPNet ===\n");
  printf("NumNets : %ld\n", netPtrs_.size());
  for(const auto& net : netPtrs_)
    printNetInfo(net);
}

void
SkyPlaceDB::printCellInfo(const GPCell* cell) const
{
  std::cout << "GPCell ID: " << cell->id() << std::endl;
  printf("(%f, %f) - (%f, %f)\n", 
      cell->lx(), cell->ly(), cell->ux(), cell->uy() );

  printf("GPPin Coordinate\n");
  for(const auto& pin: cell->pins())
    printf("ID: %d (%f, %f) GPNet: %d\n", pin->id(), pin->cx(), pin->cy(), pin->net()->id());
}

void
SkyPlaceDB::printNetInfo(const GPNet* net) const
{
  assert(net != nullptr);
  printf("NetID : %d ", net->id());

  if(!net->pins().empty())
    printf("(%f, %f) - (%f, %f)", net->lx(), net->ly(), net->ux(), net->uy() );
  printf("\n");

  for(const auto& pin: net->pins())
  {
    assert(pin != nullptr);
    assert(pin->cell() != nullptr);
    printf("PinID: %d (%f, %f) GPCell: %d\n", pin->id(), pin->cx(), pin->cy(), pin->cell()->id());
  }
}

} // namespace skyplace
