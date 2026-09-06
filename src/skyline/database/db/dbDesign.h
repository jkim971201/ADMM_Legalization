#ifndef DB_DESIGN_H
#define DB_DESIGN_H

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#include "def/defrReader.hpp"
#include "def/defiAlias.hpp"

namespace db
{

class dbTypes;
class dbTech;
class dbDie;
class dbRow;
class dbInst;
class dbNet;
class dbBTerm;
class dbITerm;
class dbViaMaster;
class dbNonDefaultRule;
class dbBlockage;
class dbTrackGrid;
class dbLayer;
class dbRegion;
class dbGroup;

class dbDesign
{
  public:

    dbDesign(const std::shared_ptr<dbTypes> types,
             const std::shared_ptr<dbTech>  tech);
    ~dbDesign();

    void writeDef       (const char* path = "") const;
    void writeBookShelf (const char* path, int dbu_bookshelf = 1) const;
    void finish(); // post-processing after read .def

    // StartCbk
    void componentStartCbk(int num_inst);
    void pinStartCbk(int num_pin);

    // EndCbk
    void componentEndCbk();
    void pinEndCbk();

    // Setters
    void setName(const char* name) { name_ = std::string(name); }
    void setDbu(int dbu);
    void setDivider(const char div);
    void setDie(const defiBox* box);
    void setDieFromBookShelf(int lx, int ly, int ux, int uy);

    // Getters (By Name or sth)
    // FIXME : All those std::string should be changed to std::string_view
    dbInst*           getInstByName           (const std::string& name);
    dbNet*            getNetByName            (const std::string& name);
    dbBTerm*          getBTermByName          (const std::string& name);
    dbViaMaster*      getViaMasterByName      (std::string_view name);
    dbGroup*          getGroupByName          (std::string_view name);
    dbRegion*         getRegionByName         (std::string_view name);
    dbNonDefaultRule* getNonDefaultRuleByName (const std::string& name);
    dbTrackGrid*      getTrackGridByLayer     (dbLayer* layer);

    char divider() const { return divider_; }

    // TODO: This is definitely not the best to way to describe a "core region".
    int coreLx() const { return coreLx_; }
    int coreLy() const { return coreLy_; }
    int coreUx() const { return coreUx_; }
    int coreUy() const { return coreUy_; }

    void setCoreLx(int lx) { coreLx_ = lx; }
    void setCoreLy(int ly) { coreLy_ = ly; }
    void setCoreUx(int ux) { coreUx_ = ux; }
    void setCoreUy(int uy) { coreUy_ = uy; }

    // Row
    void addNewRow  (const defiRow* row);

    // Inst
    void addNewInst (const defiComponent* comp, const std::string& name);
    void fillInst   (const defiComponent* comp, dbInst* inst);

    // BTerm (IO)
    void addNewIO   (const defiPin* pin, const std::string& name);

    // Net
    dbNet* getNewNet(const std::string& name);
    void   fillNet  (const defiNet* defNet, dbNet* net);

    // Special Net
    dbNet* getNewSNet(const std::string& name);

    // Generated Via
    void addNewViaMaster(const defiVia* via);

    // NonDefaultRule
    void addNewNonDefaultRule(const defiNonDefault* ndr);
 
    // Blockage
    void addNewBlockage(const defiBlockage* blk);

    // Track
    void addNewTrack(const defiTrack* tr);

    // Region
    void addNewRegion(const defiRegion* rg);

    // Group 
    void createNewGroup(std::string_view group_name);
    void fillGroup     (const defiGroup* def_group);
    void addGroupMember(std::string_view member_name);

    // Getters
          dbDie* getDie()       { return die_;  }
    const dbDie* getDie() const { return die_;  }

          std::vector<dbRow*>& getRows()       { return row_ptrs_; }
    const std::vector<dbRow*>& getRows() const { return row_ptrs_; }

          std::vector<dbInst*>& getInsts()       { return inst_ptrs_; }
    const std::vector<dbInst*>& getInsts() const { return inst_ptrs_; }

          std::vector<dbITerm*>& getITerms()       { return iterms_; }
    const std::vector<dbITerm*>& getITerms() const { return iterms_; }

          std::vector<dbBTerm*>& getBTerms()       { return bterms_; }
    const std::vector<dbBTerm*>& getBTerms() const { return bterms_; }

          std::vector<dbNet*>& getNets()       { return nets_; }
    const std::vector<dbNet*>& getNets() const { return nets_; }

          std::vector<dbNet*>& getSNets()       { return snets_; }
    const std::vector<dbNet*>& getSNets() const { return snets_; }

          std::vector<dbViaMaster*>& getViaMasters()       { return vias_; }
    const std::vector<dbViaMaster*>& getViaMasters() const { return vias_; }

          std::vector<dbNonDefaultRule*>& getNonDefaultRules()       { return ndrs_; }
    const std::vector<dbNonDefaultRule*>& getNonDefaultRules() const { return ndrs_; }

          std::vector<dbBlockage*>& getBlockages()       { return blockages_; }
    const std::vector<dbBlockage*>& getBlockages() const { return blockages_; }

          std::vector<dbTrackGrid*>& getTrackGrids()       { return trackGrids_; }
    const std::vector<dbTrackGrid*>& getTrackGrids() const { return trackGrids_; }

          std::vector<dbRegion*>& getRegions()       { return regions_; }
    const std::vector<dbRegion*>& getRegions() const { return regions_; }

          std::vector<dbGroup*>& getGroups()       { return groups_; }
    const std::vector<dbGroup*>& getGroups() const { return groups_; }

    // Returns the design name
    const std::string& name() const { return name_; } 

  private:

    int counter_db_inst_;

    int coreLx_;
    int coreLy_;
    int coreUx_;
    int coreUy_;

    std::shared_ptr<dbTech>  tech_;
    std::shared_ptr<dbTypes> types_;

    char divider_;

    std::string name_;

    // Pointers
    dbDie* die_;
    std::vector<dbRow*>       row_ptrs_;
    std::vector<dbInst*>      inst_ptrs_;
    std::vector<dbNet*>       nets_;
    std::vector<dbNet*>       snets_;
    std::vector<dbBTerm*>     bterms_;
    std::vector<dbITerm*>     iterms_;
    std::vector<dbViaMaster*> vias_;
    std::vector<dbNonDefaultRule*> ndrs_;
    std::vector<dbBlockage*>  blockages_;
    std::vector<dbTrackGrid*> trackGrids_;
    std::vector<dbRegion*>    regions_;
    std::vector<dbGroup*>     groups_;

    // Instances
    std::vector<dbRow>        rows_;
    std::vector<dbInst>       insts_;

    std::unordered_map<std::string_view, dbGroup*> pattern2dbGroup_;
    std::unordered_map<dbLayer*, dbTrackGrid*> layer2TrackGrid_;

    std::unordered_map<std::string, dbInst*>           str2dbInst_;
    std::unordered_map<std::string, dbNet*>            str2dbNet_;
    std::unordered_map<std::string, dbBTerm*>          str2dbBTerm_;
    std::unordered_map<std::string, dbViaMaster*>      str2dbViaMaster_;
    std::unordered_map<std::string, dbRegion*>         str2dbRegion_;
    std::unordered_map<std::string, dbGroup*>          str2dbGroup_;
    std::unordered_map<std::string, dbNonDefaultRule*> str2dbNonDefaultRule_;
};

}

#endif
