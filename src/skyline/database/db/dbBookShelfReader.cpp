#include "dbBookShelfReader.h"
#include "dbMacro.h"
#include "dbLayer.h"
#include "dbRect.h"
#include "dbDie.h"
#include "dbRow.h"
#include "dbNet.h"
#include "dbInst.h"
#include "dbBTerm.h"
#include "dbITerm.h"

#include <set>
#include <unordered_map>
#include <cassert>

namespace db
{

dbBookShelfReader::dbBookShelfReader(std::shared_ptr<dbTypes>  types,
                                     std::shared_ptr<dbTech>   tech,
                                     std::shared_ptr<dbDesign> design)
  : types_        (types),
    tech_         (tech),
    design_       (design),
    dbuBookShelf_ (1)
{
  bsParser_ = std::make_unique<BookShelfParser>();
}

void
dbBookShelfReader::readFile(const std::string& filename)
{
  bsParser_->parse(filename);
  //bsParser_->printInfo();
  convert2db();
}

void
dbBookShelfReader::convert2db()
{
  // Note that BookShelf format does not have any
  // technology information.
  // Though some tech-related objects (e.g. dbMacro)
  // will be used in this function, but they are
  // just dummy objects to smoothly
  // convert bookshelf format to LEF/DEF-based database.
  const auto bookshelfDB = bsParser_->getDB();
  design_->setName(bsParser_->getBenchName()); 

  // DBU should be 1 in the Bookshelf Flow
  tech_->setDbu(1);

  // IO pins of MMS benchmark are duplicate in .nets file.
  // A hash table is necessary to avoid multiple dbBTerm for
  // an identical IO pin.
  // (NOTE 25/06/27 Multiple dbBTerm is more natural)
  std::set<BsCell*> bs_cell_to_be_bterm;
  std::unordered_map<BsCell*, dbInst*> bsCell2dbInst;

  // Bookshelf basically contains micron unit
  // so we need a scaling factor to
  // convert to integer data of dbInst, dbDie, ...
  // This is just a magic number.
  // (but this has to be an even number because
  // numbers of bookshelf files are normally
  // multiplies of 0.5)
  constexpr int dbuBookShelf = 2;
  assert(dbuBookShelf >= 2);
  dbuBookShelf_ = dbuBookShelf;

  auto convert2dbDie = [&] (BsDie* bsDie)
  {
    int die_lx = dbuBookShelf * bsDie->lx();
    int die_ly = dbuBookShelf * bsDie->ly();
    int die_ux = dbuBookShelf * bsDie->ux();
    int die_uy = dbuBookShelf * bsDie->uy();
    design_->setDieFromBookShelf(die_lx, die_ly, die_ux, die_uy);

    // Assume BookShelf Format always have same core size with die
    design_->setCoreLx(die_lx);
    design_->setCoreLy(die_ly);
    design_->setCoreUx(die_ux);
    design_->setCoreUy(die_uy);
  };

  // This will be used to determine Macro
  const int raw_row_height = bookshelfDB->rowHeight(); 
  
  int numRow = 0;
  auto convert2dbRow = [&] (BsRow* bsRow)
  {
    dbRow* newRow = new dbRow;
    const std::string row_name = "BookShelfRow" + std::to_string(numRow++);
    newRow->setName(row_name);
    newRow->setOrigX(dbuBookShelf * bsRow->lx());
    newRow->setOrigY(dbuBookShelf * bsRow->ly());
    newRow->setSiteSize(dbuBookShelf * bsRow->siteWidth(),  
                        dbuBookShelf * bsRow->dy());
    newRow->setStepX(dbuBookShelf * bsRow->siteSpacing());
    newRow->setStepY(dbuBookShelf * bsRow->dy());
    newRow->setNumSiteX(bsRow->numSites());
    newRow->setNumSiteY(1);
    return newRow;
  };

  // Since there are not MACRO information like LEF in Bookshelf,
  // we distinguish each cell by its size.
  // If two instances has width and height, then we assume
  // two instances have same Master Cell (MACRO).
  std::unordered_map<std::size_t, dbMacro*> size2macro;

  auto convert2dbInst = [&] (BsCell* bsCell)
  {
    dbInst* newInst = new dbInst;
    newInst->setName(bsCell->name());
    int macroSizeX = dbuBookShelf * bsCell->dx();
    int macroSizeY = dbuBookShelf * bsCell->dy();
    
    std::string macro_size_str 
      = std::to_string(macroSizeX) + "_" + std::to_string(macroSizeY);
    std::size_t hashing_value = std::hash<std::string>()(macro_size_str);

    dbMacro* macroPtr;
    auto itr = size2macro.find(hashing_value);
    if(itr == size2macro.end())
    {
      const std::string macro_name 
        = "BookShelfMacro" + std::to_string(size2macro.size());
      macroPtr = tech_->getNewMacro(macro_name.c_str());
      macroPtr->setSizeX(macroSizeX);
      macroPtr->setSizeY(macroSizeY);
      size2macro[hashing_value] = macroPtr;
    }
    else
      macroPtr = itr->second;
  
    // Use raw value
    if(bsCell->dy() > raw_row_height)
      macroPtr->setMacroClass(MacroClass::BLOCK);

    newInst->setMacro(macroPtr);
    newInst->setLocation(dbuBookShelf * bsCell->lx(), 
                         dbuBookShelf * bsCell->ly());

    if(bsCell->isFixed() == true)
      newInst->setStatus(Status::FIXED);

    return newInst;
  };

  // dummy layer to make dbRect
  dbLayer* dummyLayer = new dbLayer; 
  int num_bterm_made = 0;
  auto convert2dbBTerm = [&] (BsPin* bs_pin_ptr)
  {
    dbBTerm* newIO = new dbBTerm;
    BsCell* bs_cell_ptr = bs_pin_ptr->cell();

    std::string io_name = bs_cell_ptr->name() + "_IO" + std::to_string(num_bterm_made++);
    newIO->setName(io_name);

    dbBTermPort* newBTermPort = new dbBTermPort;
    newBTermPort->setOrigX(0);
    newBTermPort->setOrigY(0);  

    // Default orient is N
    const float bs_cell_cx = bs_cell_ptr->cx();
    const float bs_cell_cy = bs_cell_ptr->cy();
    const float bs_pin_offset_x = bs_pin_ptr->offsetX();
    const float bs_pin_offset_y = bs_pin_ptr->offsetY();

    const float new_pin_cx = bs_cell_cx + bs_pin_offset_x;
    const float new_pin_cy = bs_cell_cy + bs_pin_offset_y;

    newBTermPort->setLx(dbuBookShelf * new_pin_cx);
    newBTermPort->setLy(dbuBookShelf * new_pin_cy);

    newBTermPort->setUx(dbuBookShelf * new_pin_cx);
    newBTermPort->setUy(dbuBookShelf * new_pin_cy);
    newBTermPort->setLayer( dummyLayer );
    newIO->addPort( newBTermPort );

    return newIO;
  };

  std::unordered_map<dbMacro*, std::unordered_map<std::size_t, dbMTerm*>> mterm_table;

  auto& dbITermVector = design_->getITerms();
  auto& dbBTermVector = design_->getBTerms();

  for(auto& [size_hash, macro_ptr]: size2macro)
  {
    mterm_table.insert(
        std::make_pair(macro_ptr, std::unordered_map<std::size_t, dbMTerm*>()));
  }

  auto convert2dbNet = [&] (BsNet* bsNet)
  {
    dbNet* newNet = new dbNet;
    newNet->setName(bsNet->name());

    for(auto bsPinPtr : bsNet->pins())
    {
      BsCell* bsCellPtr = bsPinPtr->cell();

      bool is_io_pin = bs_cell_to_be_bterm.count(bsCellPtr) > 0 ? true : false;
      if(is_io_pin == false)
      {
        dbITerm* newITerm = new dbITerm;
        newITerm->setNet(newNet);

        // Set dbInst
        auto inst_itr = bsCell2dbInst.find(bsCellPtr);
        dbInst* inst_ptr;
        if(inst_itr == bsCell2dbInst.end())
        {
          printf("Error while converting bookshelf to db...\n");
          exit(1);
        }
        else
        {
          inst_ptr = inst_itr->second;
          const std::string iterm_name 
            = bsCellPtr->name() + "_" + std::to_string(inst_ptr->getITerms().size());
          newITerm->setName(iterm_name);
          newITerm->setInst(inst_ptr);
        }

        // Set dbMTerm
        dbMacro* macro_ptr = inst_ptr->macro();

        auto& offset2dbMTerm = mterm_table[macro_ptr];
        int offsetX = dbuBookShelf * bsPinPtr->offsetX();
        int offsetY = dbuBookShelf * bsPinPtr->offsetY();

        std::string iterm_str
          = std::to_string(offsetX) + "_" 
          + std::to_string(offsetY) + "_"
          + std::to_string(inst_ptr->getITerms().size());
        std::size_t hashing_value_iterm = std::hash<std::string>()(iterm_str);

        dbMTerm* mterm;
        auto mterm_itr = offset2dbMTerm.find(hashing_value_iterm);
        if(mterm_itr == offset2dbMTerm.end())
        {
          mterm = new dbMTerm;
          dbMTermPort* newPort = new dbMTermPort;
          newPort->setLayer(dummyLayer); // this is just a dummy layer
          newPort->addPoint(offsetX, offsetY);
          // POLYGON implicitly finish with the starting point.
          mterm->addPort( newPort );
          mterm->setBoundary();
          const std::string mterm_name 
           = macro_ptr->name() + "_" + iterm_str;
          mterm->setName(mterm_name);
          macro_ptr->addMTerm(mterm);
          offset2dbMTerm[hashing_value_iterm] = mterm;
        }
        else
          mterm = mterm_itr->second;

        newITerm->setMTerm(mterm);

        inst_ptr->addITerm(newITerm); 
        // addIterm checks if dbMTerm exists in dbITerm.
        // so this must be called after assign dbMTerm of this ITerm.

        // Finish ITerm
        newNet->addITerm(newITerm);
        dbITermVector.push_back(newITerm); 
        // new dbITerm should be added to dbDatabse
      }
      else
      {
        auto new_bterm_ptr = convert2dbBTerm(bsPinPtr);
        dbBTermVector.push_back(new_bterm_ptr);
        newNet->addBTerm(new_bterm_ptr);
        new_bterm_ptr->setNet(newNet);
      }
    }

    return newNet;
  };

  // Die
  convert2dbDie(bookshelfDB->getDie());

  // Row
  auto& dbRowVector = design_->getRows();
  auto& bsRowVector = bookshelfDB->rowVector();
  for(auto bsRowPtr : bsRowVector)
    dbRowVector.push_back(convert2dbRow(bsRowPtr));

  // Inst && BTerm (Bookshelf describes both IOs and instances in the same way)
  auto& dbInstVector = design_->getInsts();
  auto& bsCellVector = bookshelfDB->cellVector();
  for(auto bsCellPtr : bsCellVector)
  {
    // Bookshelf format does not have IOs explicitly.
    // we have to detect them by context.
    // (I think "terminal" keyword does not mean that it is an IO pin.)
    if(bsCellPtr->dx() == 0 || bsCellPtr->dy() == 0 || 
      (bsCellPtr->isFixed() && bsParser_->isOutsideDie(bsCellPtr)) )
    {
      // Record BookShelfCells, which is going to be dbBTerm (IO)
      bs_cell_to_be_bterm.insert(bsCellPtr);
    }

    auto newInst = convert2dbInst(bsCellPtr);
    dbInstVector.push_back(newInst);
    bsCell2dbInst[bsCellPtr] = newInst;
  }

  // Net
  auto& bsNetVector = bookshelfDB->netVector();
  auto& dbNetVector = design_->getNets();
  for(auto bsNetPtr : bsNetVector)
  {
    auto newNet = convert2dbNet(bsNetPtr);
    dbNetVector.push_back(newNet);
  }

  // Though bookshelf format does not have 
  // technology information,
  // we should add dummy layer information for 
  // graphic interface.
  tech_->addNewLayerBookShelf(dummyLayer);

  printf("  Finish DB converting\n");
}

}
