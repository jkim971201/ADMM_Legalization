#ifndef SKYPLACE_DB_H
#define SKYPLACE_DB_H

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#define SQRT2 1.414213562373095048801L
#define DENSITY_SCALE SQRT2

namespace db 
{
  class dbDatabase;
  class dbTech;
  class dbDesign;
  class dbTypes;
  class dbDie;
  class dbInst;
  class dbNet;
  class dbITerm;
  class dbBTerm;
}

namespace skyplace 
{

using namespace db;
using OverlapBins = std::pair<std::pair<int, int>, std::pair<int, int>>;

class GPPin;
class GPNet;
class GPCell;
class GPBin;
class GPDie;

class SkyPlaceDB 
{
  public: 

    SkyPlaceDB(); // Constructor (Just for Initialization)

    // Important APIs
    void reset();
    void exportDB(std::shared_ptr<dbDatabase> _dbDatabase); 
    void init(std::shared_ptr<dbDatabase> db); // Initialization to run placement 
    void setTargetDensity (float density) { targetDensity_ = density; }
    void setNumBinX       (int   numBinX) { numBinX_       = numBinX; }     
    void setNumBinY       (int   numBinY) { numBinY_       = numBinY; }     

    void updateMacroDensityWeight(float macroWeight);

    // Used by Initial Placer
    void moveCellInsideLayout(GPCell* cell); 
    void setNumCluster(int numCluster) { numCluster_  = numCluster; }

    // Getters
    const std::string&   designName() const { return designName_;  }
    const std::string&   designDir () const { return  designDir_;  }
    const std::vector<GPCell*>& cells() const { return   cellPtrs_;  }
    const std::vector<GPNet*>&   nets() const { return    netPtrs_;  }
    const std::vector<GPPin*>&   pins() const { return    pinPtrs_;  }
    const std::vector<GPBin*>&   bins() const { return    binPtrs_;  }
    const std::vector<GPPin*>&   getIOPins() const { return io_pins_; }

    const std::vector<GPCell*>&   fixedCells() const { return fixedPtrs_;   }
    const std::vector<GPCell*>& movableCells() const { return movablePtrs_; }

    std::shared_ptr<GPDie> die() const { return die_ptr_; }

    int numFixed   () const { return numFixed_;        }
    int numMovable () const { return numMovable_;      }
    int numFiller  () const { return numFiller_;       }
    int numNet     () const { return netPtrs_.size();  }
    int numPin     () const { return pinPtrs_.size();  }
    int numIO      () const { return numIO_;           }
    int numCluster () const { return numCluster_;      }
    int numMacro   () const { return numMovableMacro_; }
    int getDbu     () const { return dbu_;             }
    int getBookshelfDbu() const { return bookshelf_dbu_; }

    float util                 () const { return util_;                 }
    float density              () const { return density_;              }
    float targetDensity        () const { return targetDensity_;        }
    float sumFixedArea         () const { return sumFixedArea_;         }
    float sumMovableArea       () const { return sumMovableArea_;       }
    float sumScaledMovableArea () const { return sumScaledMovableArea_; }
    // sumMovableArea does not include fillerArea
    // sumScaledMovableArea is used for createFillers()
    // sumScaledMovableArea = sumStdCells + sumMacroArea * targetDensity

    float getHPWL () const { return hpwl_;    }
    int   numBinX () const { return numBinX_; } // ex) 512x256 Grid -> numBinX = 512
    int   numBinY () const { return numBinY_; } // ex) 512x256 Grid -> numBinY = 256
    float binX    () const { return binX_;    }
    float binY    () const { return binY_;    }

    void  updateHpwl();
    void  updatePinBound(); // For B2B Model (CG-based Initialization)
  
    // To plot density gradient arrows, these will be delivered to Painter
    std::vector<float>& densityGradX() { return densityGradX_; }
    std::vector<float>& densityGradY() { return densityGradY_; }

    // For logging
    void printInfo() const;

    // For Debugging
    void debugCell() const;
    void debugNet () const;

  private:

    std::string designName_;
    std::string designDir_;

    float targetDensity_;

    int dbu_;           // This is from dbDatabase->dbTech->getDbu()
    int bookshelf_dbu_; // This is from dbDatabase->getBookshelfDbu()
    bool is_bookshelf_;

    // Number of Objects
    int numStdCells_;
    int numMacro_;
    int numFixed_;
    int numMovable_;
    int numMovableMacro_;
    int numFiller_;
    int numBinX_;
    int numBinY_;
    int numIO_;
    int numCluster_;

    // Sub-Routines of init                         
    void validateTargetDensity();
    void importDB(std::shared_ptr<dbDatabase> _dbDatabase); // Step #1
    void createBins();                                      // Step #2
    void createFillers();                                   // Step #3
    void updateDensitySize();                               // Step #4
    void updateFixedOverlapArea();                          // Step #5
    int  numInitStep_;                                      // For logging

    // GPBin-related Methods
    OverlapBins findBin                  (const GPCell* cell);
    OverlapBins findBinWithDensitySize   (const GPCell* cell);

    bool isOutsideDie                    (const GPCell* cell);
    bool isOverlap                       (const GPCell* cell, const GPBin* bin);

    float getXCoordiInsideLayout         (const GPCell* cell);
    float getYCoordiInsideLayout         (const GPCell* cell);

    float getXCoordiInsideLayout         (const GPCell* cell, float x);
    float getYCoordiInsideLayout         (const GPCell* cell, float y);

    float getXDensityCoordiInsideLayout  (const GPCell* cell);
    float getYDensityCoordiInsideLayout  (const GPCell* cell);

    // GPBin-related
    float sumFixedArea_;
    float sumMovableArea_;
    float sumMovableStdArea_;
    float sumMovableMacroArea_;
    float sumScaledMovableArea_;

    float sumTotalInstArea_;
    float dieArea_;
    float density_;
    float util_;

    float avgCellArea_;
    float avgStdCellArea_;
    float idealBinAreaForAvgCellArea_;
    float idealBinAreaForAvgStdCellArea_;
    float idealBinCountForAvgCellArea_;
    float idealBinCountForAvgStdCellArea_;
    float idealNumBinForAvgCellArea_;
    float idealNumBinForAvgStdCellArea_;

    int   idealNumBin_;
    float binX_;
    float binY_;

    float avgFillerArea_;
    float whiteSpaceArea_;
    float fillerArea_;

    float fillerWidth_;
    float fillerHeight_;

    float hpwl_;

    // Database
    std::vector<GPCell*> cellPtrs_; 
    std::vector<GPCell>  cellInsts_; 

    std::vector<GPCell*> movableMacroPtrs_; 
    std::vector<GPCell*> fixedPtrs_; 
    std::vector<GPCell*> movablePtrs_; 

    std::vector<GPNet*>  netPtrs_; 
    std::vector<GPNet>   netInsts_; 

    std::vector<GPPin*>  io_pins_;
    std::vector<GPPin*>  pinPtrs_; 
    std::vector<GPPin>   pinInsts_; 

    std::vector<GPBin*>  binPtrs_;
    std::vector<GPBin>   binInsts_;

    std::unordered_map<dbInst*,  GPCell*> dbInst2Cell_;

    std::shared_ptr<GPDie> die_ptr_;

    // To plot density gradient arrows...
    std::vector<float> densityGradX_;
    std::vector<float> densityGradY_;

    // For Debug
    void printBinDensity();
    void printBinPotential();
    void printBinElectroFieldXY();
    void printFixedOverlapArea();

    void printCellInfo(const GPCell* cell) const;
    void printNetInfo(const GPNet* net) const;
};

} // namespace skyplace 

#endif
