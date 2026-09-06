#include <iostream>
#include <cassert> 
#include <cstring>
#include <filesystem>
#include <stdio.h>
#include <string.h>

// For C++-Style Parsing
// 2025/03/04
#include <vector>
#include <fstream>
#include <string>
#include <string_view>

#include "BookShelfParser.h"

#define BUF_SIZE 512

namespace bookshelf
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

inline char* getNextWord(const char* delimiter)
{
  // strtok() has "static char *olds" inside
  // if strtok() gets nullptr, then it starts from the "olds"
  // so you can get the next word
  // you must use strtok(buf, ~) to capture the first word
  // otherwise, you will get segmentation fault
  // (this is implemented in goNextLine)
  return strtok(nullptr, delimiter);
}

// Go Next Line and Get First Word
inline char* goNextLine(char* buf, const char* delimiter, FILE* fp)
{
  fgets(buf, BUF_SIZE-1, fp);
  return strtok(buf, delimiter);
}

inline void extractDir(char* file, char* dir)
{
  int nextToLastSlash = 0;
  int len = strlen(file);

  strcpy(dir, file);

  for(int i = 0; i < len; i++) 
  {
    if(file[i] == '/') 
      nextToLastSlash = i + 1;
  }

  dir[nextToLastSlash] = '\0';
  // \0 == nullptr
}

// What a stupid way...
// But this seems best... at least to me...
inline void catDirName(char* dir, char* name)
{
  char temp[MAX_FILE_NAME];
  strcpy(temp, dir);
  strcat(temp, name);
  strcpy(name, temp);
}

inline void getOnlyName(const char* file, char* name)
{
  int lastDot = 0;
  int nextToLastSlash = 0;
  int len = strlen(file);

  for(int i = 0; i < len; i++) 
  {
    if(file[i] == '.') 
      lastDot = i;
    if(file[i] == '/') 
      nextToLastSlash = i + 1;
  }

  int nameLength = lastDot - nextToLastSlash;
  for(int i = 0; i < nameLength; i++)
    name[i] = file[i + nextToLastSlash];

  name[nameLength] = '\0';
}

inline void getSuffix(const char* token, char* sfx)
{
  int lastDot = 0;
  int len = strlen(token);

  for(int i = 0; i < len; i++) 
  {
    if(token[i] == '.') 
      lastDot = i;
  }
  strcpy(sfx, &token[lastDot + 1]);
}

// BsCell //
BsCell::BsCell()
{
  lx_ = ly_ = ux_ = uy_ = 0;
  id_ = 0;
  orient_ = 'N';
  isFixed_ = isFixedNI_ = false;
}

BsCell::BsCell(std::string& name,
               int  id,
               int  width, 
               int  height, 
               bool isTerminal, 
               bool isTerminalNI) 
    : BsCell()
{
  name_ = name;

  id_ = id;
  dx_ = width;
  dy_ = height;
  isTerminal_   = isTerminal;
  isTerminalNI_ = isTerminalNI;
}

// BsRow //
BsRow::BsRow(int idx,
             int ly, 
             int rowHeight, 
             int siteWidth, 
             int siteSpacing,
             int offsetX,
             int numSites)
{
  idx_          = idx;

  // Explicit Values from .scl file
  ly_           = ly;
  rowHeight_    = rowHeight;
  siteWidth_    = siteWidth;
  siteSpacing_  = siteSpacing;
  offsetX_      = offsetX;
  numSites_     = numSites;
  siteOrient_   = true;
  siteSymmetry_ = true;

  // Implicit Value
  rowWidth_ = numSites_ * siteWidth_;
}

// BsPin //
BsPin::BsPin(BsCell* cell, int netID, 
             double offsetX, double offsetY,
             char IO)
{
  cell_  = cell;
  net_   = nullptr;

  netID_ = netID;

  offsetX_ = offsetX;
  offsetY_ = offsetY;

  io_ = IO;
}

// BookShelfDB //
BookShelfDB::BookShelfDB(int numNodes)
{
  numRows_    = 0;
  numFixed_   = 0;
  numFixedNI_ = 0;
  numMovable_ = 0;
  numInst_    = 0;
  rowHeight_  = 0;

  numBsCells_ = numNodes;
  cellPtrs_.reserve(numNodes);
  cellInsts_.reserve(numNodes);
}

void
BookShelfDB::makeBsCell(std::string& name, 
                        int  width, 
                        int  height, 
                        bool isTerminal, 
                        bool isTerminalNI)
{
  BsCell oneBsCell(name, numInst_, 
                   width, height, 
                   isTerminal, isTerminalNI);

  numInst_++;

  // TODO: Correct this
  if(isTerminal)        numFixed_++;
  else if(isTerminalNI) numFixedNI_++;
  else                  numMovable_++;

  cellInsts_.push_back(oneBsCell);
  // cellPtrs will be filled by buildBsCellMap()
}

void
BookShelfDB::buildBsCellMap()
{
  //printf("  Building Node Map\n");
  for(BsCell& c : cellInsts_)
  {
    cellPtrs_.push_back(&c);
    cellMap_.emplace(c.name(), &c);
  }
}

void
BookShelfDB::makeBsRow(int idx,
                       int ly, 
                       int rowHeight, 
                       int siteWidth, 
                       int siteSpacing,
                       int offsetX,
                       int numSites)
{
  BsRow oneBsRow(idx, ly, 
                 rowHeight, 
                 siteWidth, 
                 siteSpacing, 
                 offsetX, 
                 numSites);
  rowInsts_.push_back(oneBsRow);
  // rowPtrs will be filled by buildBsRowMap()
}

void
BookShelfDB::makeBsNet(int netID, const std::string& name)
{
  BsNet oneBsNet(netID);
  oneBsNet.setName(name);
  netInsts_.push_back(oneBsNet);
  // netPtrs will be filled by finishPinsAndNets()
}

void
BookShelfDB::makeBsPin(BsCell* cell, int netID, 
                       double offsetX, double offsetY,
                       char IO)
{
  double cell_half_w = static_cast<double>(cell->dx()) / 2;
  double cell_half_h = static_cast<double>(cell->dy()) / 2;

//  if(std::abs(offsetX) > cell_half_w)
//  {
//    printf("  [Warning] Pin of Cell %s is out of Cell Boundary.\n", 
//              cell->name().c_str());
//    printf("  [Warning] Pin OffsetX %.1f is larger than Cell Half Width %.1f \n",
//              offsetX, cell_half_w);
//  }
//
//  if(std::abs(offsetY) > cell_half_h)
//  {
//    printf("  [Warning] Pin of Cell %s is out of Cell Boundary.\n", 
//              cell->name().c_str());
//    printf("  [Warning] Pin OffsetY %.1f is larger than Cell Half Height %.1f \n",
//              offsetY, cell_half_h);
//  }

  BsPin oneBsPin(cell, netID, offsetX, offsetY, IO);

  pinInsts_.push_back(oneBsPin);
}

void
BookShelfDB::finishPinsAndNets()
{
  //printf("Building Net maps.\n");

  for(auto& net : netInsts_)
  {
    netPtrs_.push_back(&net); 
    netMap_[net.id()] = &net; 
  }

  //printf("  Adding pins to cells and nets.\n");

  for(auto& pin : pinInsts_)
  {
    pinPtrs_.push_back(&pin);
    pin.cell()->addNewPin(&pin);
    if(!pin.net()) 
      pin.setNet(getBsNetByID(pin.netID()));
    else
    {
      printf("  Unknown Error...\n");
      exit(0);
    }

    if(pin.net()) 
      pin.net()->addNewPin(&pin);
    else
    {
      printf("  Unknown Error...\n");
      exit(0);
    }
  }
}

void
BookShelfDB::buildBsRowMap()
{
  int maxX = 0;
  int maxY = 0;

  int minX = std::numeric_limits<int>::max();
  int minY = std::numeric_limits<int>::max();

  //printf("  Building Row Map\n");
  for(BsRow& r : rowInsts_)
  {
    if(r.ux() > maxX) maxX = r.ux();
    if(r.lx() < minX) minX = r.lx();
    if(r.uy() > maxY) maxY = r.uy();
    if(r.ly() < minY) minY = r.ly();
    rowPtrs_.push_back(&r);
    rowMap_.emplace(r.id(), &r);
  }

  //for(auto& c : cellPtrs_)
  //  if(c->uy() < minY) minY = c->uy();

  bsDie_.setUxUy(maxX, maxY);
  bsDie_.setLxLy(minX, minY);
  bsDiePtr_ = &bsDie_;

  //printf("Creating a Die (%d, %d) - (%d, %d) \n", maxX, maxY, minX, minY);
  numRows_ = rowPtrs_.size();
}

void
BookShelfDB::verifyMap()
{
  std::cout << "Start Verifying Map Vector" << std::endl;
  for(auto& kv : cellMap_)
  {
    std::cout << "key name: " << kv.first << std::endl;
    std::cout << "ptr width: " << kv.second->dx() << std::endl;
  }
}

void
BookShelfDB::verifyVec()
{
  std::cout << "Start Verifying Instance Vector" << std::endl;
  for(auto& c : cellInsts_)
  {
    std::cout << "cell name: " << c.name() << std::endl;
    std::cout << "cell width: " << c.dx() << std::endl;
  }
}

void
BookShelfDB::verifyPtrVec()
{
  std::cout << "Start Verifying Pointer Vector" << std::endl;
  for(int i = 0; i < cellPtrs_.size(); i++)
  {
    std::cout << "cell name: " << cellPtrs_[i]->name() << std::endl;
    std::cout << "cell width: " << cellPtrs_[i]->dx() << std::endl;
  }
}

BookShelfParser::BookShelfParser()
{
  numFixed_     = 0;
  numFixedNI_   = 0;
  numMovable_   = 0;
  maxRowHeight_ = 0;
  bookShelfDB_  = nullptr;
}

void
BookShelfParser::init(const char* aux_name)
{
  char suf[4];
  getSuffix(aux_name, suf);

  getOnlyName(aux_name, benchName_);

  if(!strcmp(suf, "aux"))
  {
    strcpy(aux_, aux_name);
    extractDir(aux_, dir_);
  }
  else
  {
    printf("Make sure you give .aux file\n");
    exit(0);
  }
}

void
BookShelfParser::parse(const std::filesystem::path& aux_name)
{
  if(!std::filesystem::exists(aux_name))
  {
    printf("aux file %s not exists!\n", aux_name.c_str());
    exit(1);
  }

  // filetype check will be done in BookShelfParser::init
  std::string filename = std::string(aux_name);
  init(filename.c_str());

  // Start from .aux file
  read_aux();
  read_nodes();
  read_pl();
  read_scl();
  read_nets();

  // printf(" Parsing is finished successfully!\n");
}

void
BookShelfParser::read_aux()
{
  //printf("  Reading %s...\n", aux_);
  FILE *fp = fopen(aux_, "r");
  char *token = nullptr;
  char buf[BUF_SIZE-1];

  char sfx[6]; // the longest is "nodes" (5 letters)

  if(fp == nullptr)
  {
    printf(" Failed to open %s...\n", aux_);
    exit(0);
  }
  
  token = goNextLine(buf, " :", fp);

  if(strcmp(token, "RowBasedPlacement"))
  {
    printf(" Unknown Placement Type: %s\n", token);
    exit(0);
  }

  while(true)
  {
    token = getNextWord(" :\n");
    if(!token) break;  

    getSuffix(token, sfx);

    if(!strcmp(sfx, "nodes"))
    {
     //printf(" .nodes file detected.\n");
     strcpy(nodes_, token);
     catDirName(dir_, nodes_);
    }
    else if(!strcmp(sfx, "pl"))
    {
     //printf(" .pl file detected.\n");
     strcpy(pl_, token);
     catDirName(dir_, pl_);
    }
     else if(!strcmp(sfx, "scl"))
    {
      //printf(" .scl file detected.\n");
      strcpy(scl_, token);
      catDirName(dir_, scl_);
    }
    else if(!strcmp(sfx, "nets"))
    {
      //printf(" .nets file detected.\n");
      strcpy(nets_, token);
      catDirName(dir_, nets_);
    }
    else if(!strcmp(sfx, "wts"))
    {
      printf("  .wts file is not supported. (will be ignored)\n");
    }
    else
    {
      printf("Unknown file format %s...\n", sfx);
      exit(1);
    }
  }

  fclose(fp);
}

void
BookShelfParser::read_nodes()
{
  printf("  Reading %s...\n", nodes_);

  std::ifstream input_file(nodes_);

  std::string line = std::string();
  std::string delim = ""; 

  int num_nodes = 0;
  int num_terminals = 0;

  bool find_num_nodes = false;
  bool find_num_terminals = false;

  while(std::getline(input_file, line))
  {
    auto tokens = tokenize(line, delim);
    if(tokens.empty() == true)
      continue;

    // Skip Headlines
    if(tokens.at(0) == "UCLA" || tokens.at(0) == "#")
      continue;

    if(tokens[0] == "NumNodes")
    {
      num_nodes = std::stoi(tokens.at(2));
      find_num_nodes = true;
    }
    if(tokens[0] == "NumTerminals")
    {
      num_terminals = std::stoi(tokens.at(2));
      find_num_terminals = true;
    }

    if(find_num_nodes == true && find_num_terminals == true)
      break;
  }

  if(find_num_nodes == false || find_num_terminals == false)
  {
    printf("  Failed to find NumNodes or NumTerminals...\n");
    exit(0);
  }

  bookShelfDB_ = std::make_shared<BookShelfDB>(num_nodes);
  // printf("NumNodes : %d NumTerminals : %d\n", num_nodes, num_terminals);

  std::string cell_name = "";
  int cell_width = 0;
  int cell_height = 0;
  while(std::getline(input_file, line))
  {
    bool is_terminal = false;
    bool is_terminalNI = false;
    auto tokens = tokenize(line, delim);
    if(tokens.empty() == true)
      continue;

    // Skip comments
    if(tokens.at(0) == "#")
      continue;
   
    if(tokens.size() == 3 || tokens.size() == 4)
    {
      cell_name = tokens.at(0);
      cell_width = std::stoi(tokens.at(1));
      cell_height = std::stoi(tokens.at(2));

      if(tokens.size() == 4)
      {
        if(tokens.at(3) == "terminal")
          is_terminal = true;
        else if(tokens.at(3) == "terminal_NI")
          is_terminalNI = true;
        else
        {
          printf("  Unknown Keyword %s\n", tokens.at(3).c_str());
          exit(0);
        }
      }

      // printf("CellName : %s\n", cell_name.c_str());

      bookShelfDB_->makeBsCell(cell_name, 
                               cell_width, 
                               cell_height, 
                               is_terminal, 
                               is_terminalNI);
    }
    else
    {
      // all the lines should be "cell_name cell_width cell_height"
      assert(0);
    }
  }

  bookShelfDB_->buildBsCellMap();
  // printf("NumNodes : %d CellVector : %ld\n", num_nodes, bookShelfDB_->cellVector().size());
  assert(num_nodes == bookShelfDB_->cellVector().size());
}

void
BookShelfParser::read_pl()
{
  printf("  Reading %s...\n", pl_);

  std::ifstream input_file(pl_);

  std::string line = std::string();
  std::string delim = ""; 

  int numLines = 0;

  float lx = 0;
  float ly = 0;
  
  std::string orient;

  bool isFixed   = false;
  bool isFixedNI = false;

  while(std::getline(input_file, line))
  {
    // Skip empty line
    if(line.length() <= 1)
      continue;

    // Skip Comment lines
    if(line[0] == '#')
      continue;

    // Skip UCLA Headline
    auto tokens = tokenize(line, delim);
    if(tokens[0] == "UCLA")
      continue;

    assert(tokens.size() == 5 || tokens.size() == 6);

    std::string cellName = tokens[0];

    BsCell* myBsCell = bookShelfDB_->getBsCellByName(cellName);
    assert(cellName == myBsCell->name());

    // Get X Coordinate
    lx = std::stof(tokens[1]);

    // Get Y Coordinate
    ly = std::stof(tokens[2]);

    myBsCell->setXY(lx, ly);

    assert(tokens[3] == ":");

    // Get Orient
    orient = tokens[4];

    if(orient != "N") 
      myBsCell->setOrient('N');

    // Get Move-Type
    if(tokens.size() == 6)
    {
      if(tokens[5] == "/FIXED")
      {
        myBsCell->setFixed();
        numFixed_++;
      }
      else if(tokens[5] == "/FIXED_NI")
      {
        myBsCell->setFixedNI();
        numFixedNI_++;
      }
      else
        numMovable_++;
    }

    // printf(" %s: %f %f\n", cellName.c_str(), lx, ly);
  
    numLines++;
    //if(numLines % 100000 == 0)
    //  printf("  Completed %d lines\n", numLines);
  }
}


void
BookShelfParser::read_scl()
{
  printf("  Reading %s...\n", scl_);

  std::ifstream input_file(scl_);

  std::string line = std::string();
  std::string delim = ""; 

  // Skip Headlines
  int num_rows = 0;
  int rows_read = 0;
  int max_row_height = 0;

  auto readOneRow = [&] ()
  {
    std::string line_row;

    int ly;           
    int row_height;
    int site_width;    
    int site_spacing;  
    std::string site_orient;   
    std::string site_symmetry;
    int offset_x;
    int num_sites;     

    while(std::getline(input_file, line_row))
    {
      // Skip Comment lines
      if(line_row[0] == '#')
        continue;

      auto tokens = tokenize(line_row, delim);
      if(tokens[0] == "Coordinate")
      {
        ly = std::stoi(tokens[2]);
        //printf("Ly : %d\n", ly);
      }
      else if(tokens[0] == "Height")
      {
        row_height = std::stoi(tokens[2]);
        //printf("row_height : %d\n", row_height);
        if(row_height > max_row_height)
          max_row_height = row_height;
      }
      else if(tokens[0] == "Sitewidth")
      {
        site_width = std::stoi(tokens[2]);
        //printf("site_width : %d\n", site_width);
      }
      else if(tokens[0] == "Sitespacing")
      {
        site_spacing = std::stoi(tokens[2]);
        //printf("site_spacing : %d\n", site_spacing);
      }
      else if(tokens[0] == "Siteorient")
      {
        site_orient = tokens[2];
        //printf("site_orient : %s\n", site_orient.c_str());
      }
      else if(tokens[0] == "Sitesymmetry")
      {
        site_symmetry = tokens[2];
        //printf("site_symmetry : %s\n", site_symmetry.c_str());
      }
      else if(tokens[0] == "SubrowOrigin")
      {
        offset_x = std::stoi(tokens[2]);
        //printf("offset_x : %d\n", offset_x);
        assert(tokens[3] == "NumSites");
        num_sites = std::stoi(tokens[5]);
        //printf("num_sites : %d\n", num_sites);
      }
      else if(tokens[0] == "End")
      {
        // Strictly, END is not allowed in Bookshelf format.
        break;
      }
      else
      {
        printf("  Wrong BookShelf Syntax\n");
        printf("  Current Line : %s\n", line_row.c_str());
        exit(0);
      }
    }

    bookShelfDB_->makeBsRow(rows_read, // idx
                            ly, 
                            row_height, 
                            site_width, 
                            site_spacing, 
                            offset_x, 
                            num_sites);
    rows_read++;
  };


  while(std::getline(input_file, line))
  {
    // Skip empty line
    if(line.length() <= 1)
      continue;

    // Skip Comment lines
    if(line[0] == '#')
      continue;

    // Skip UCLA Headline
    auto tokens = tokenize(line, delim);
    if(tokens[0] == "UCLA")
      continue;

    if(tokens[0] == "NumRows")
    {
      assert(tokens.size() == 3);
      num_rows = std::stoi(tokens[2]);
      // printf("NumRows : %d\n", num_rows);
    }
    else if(tokens[0] == "CoreRow")
    {
      if(tokens[1] != "Horizontal")
      {
        printf("  No support for Vertical Row!\n");
        exit(0);
      }
      readOneRow();
    }

    if(rows_read > num_rows) 
    {
      printf("  Extra Rows more than %d will be ignored...\n", num_rows);
      break;
    }
  }

  bookShelfDB_->buildBsRowMap();
  bookShelfDB_->setHeight(max_row_height);
  assert(num_rows == bookShelfDB_->rowVector().size());
  //printf("  Successfully Finished %s!\n", scl_);
}

void
BookShelfParser::read_nets()
{
  printf("  Reading %s...\n", nets_);

  FILE *fp = fopen(nets_, "r");
  char *token = nullptr;
  char buf[BUF_SIZE-1];

  if(fp == nullptr)
  {
    printf("Failed to open %s...\n", nets_);
    exit(0);
  }

  // Skip Headlines
  while(!token || !strcmp(token, "UCLA") || token[0] == '#' || std::string(token).length() <= 1)
    token = goNextLine(buf, " \t\n", fp);

  // Read NumNets (at this moment, buf == "NumNets")
  assert(!strcmp(buf, "NumNets"));
  token = getNextWord(" \t\n:");
  int numNets = atoi(token);
  token = goNextLine(buf, " \t\n", fp);

  // Read NumPins (at this moment, buf == "NumPins")
  assert(!strcmp(buf, "NumPins"));
  token = getNextWord(" \t\n:");
  int numPins = atoi(token);
  token = getNextWord(" \t\n");
  token = goNextLine(buf, " \t\n", fp);

  //printf(" Total Nets: %d\n", numNets);
  //printf(" Total Pins: %d\n", numPins);

  // Go to Next Line untill there are no blank lines anymore
  while(!token || token[0] == '#' || std::string(token).length() <= 1)
    token = goNextLine(buf, " \t\n", fp);

  int netsRead = 0;

  while(!feof(fp))
  {
    // Read a Net (at this moment, token == "NetDegree")
    assert(!strcmp(token, "NetDegree"));
    token = getNextWord(" \t\n:");
    int netDegree = atoi(token);

    token = getNextWord(" \t\n");
    std::string netName = std::string(token);

    // Some nets have NetDegree : 1 in MMS benchmarks
    // Should we ignore these netDegree 1 nets?
    bookShelfDB_->makeBsNet(netsRead, netName); // netsRead => netID

    char IO;
    double offsetX;
    double offsetY;

    int pinsRead = 0;
    while(pinsRead < netDegree)
    {
      // Get Master Cell's Name
      token = goNextLine(buf, " \t", fp);
      BsCell* cellOfThesePins = bookShelfDB_->getBsCellByName(std::string(token));

      // Get Pin IO Type
      token = getNextWord(" \t\n");
      IO = token[0];

      if(IO != 'I' && IO != 'O' && IO != 'B')
      {
        printf("  Wrong BookShelf Syntax.\n");
        printf("  Pin IO type must be one of I, O, B.\n");
        exit(0);
      }

      // Get Pin Offset X
      token = getNextWord(" \t\n:");
      if(token == nullptr) // In ICCAD 2004 Benchmarks,
        offsetX = 0.0;  // They do not provide offsets of Chip IOs
      else 
        offsetX = atof(token);

      // Get Pin Offset Y
      token = getNextWord(" \t\n:");
      if(token == nullptr) // In ICCAD 2004 Benchmarks,
        offsetY = 0.0;  // They do not provide offsets of Chip IOs
      else 
        offsetY = atof(token);

      // Ignore pins of Degree-1 net?
      bookShelfDB_->makeBsPin(cellOfThesePins, 
                              netsRead, // netID 
                              offsetX, offsetY, IO);
      pinsRead++;
    }

    netsRead++;
    //if(netsRead % 100000 == 0)
    //  printf("  Completed %d nets\n", netsRead);

    token = goNextLine(buf, " \t\n", fp);
  }

  bookShelfDB_->finishPinsAndNets();
  //printf("Successfully Finished %s!\n", nets_);
}

// Temporary...
bool 
BookShelfParser::isOutsideDie(BsCell* cell)
{
  bool isOutside = false;

  BsDie* die = bookShelfDB_->getDie();

  int dieLx = die->lx();
  int dieLy = die->ly();
  int dieUx = die->ux();
  int dieUy = die->uy();

  if(cell->ux() <= dieLx) isOutside = true;
  if(cell->lx() >= dieUx) isOutside = true;
  if(cell->uy() <= dieLy) isOutside = true;
  if(cell->ly() >= dieUy) isOutside = true;

  return isOutside;
}

void
BookShelfParser::printInfo()
{
  using namespace std;

  string name = string(benchName_);

  int numInst = numFixed_ + numMovable_;
  int numNet  = bookShelfDB_->netVector().size();
  int numPin  = bookShelfDB_->pinVector().size();
  int numRow  = bookShelfDB_->rowVector().size();

  BsDie* die = bookShelfDB_->getDie();

  int dieLx = die->lx();
  int dieLy = die->ly();
  int dieUx = die->ux();
  int dieUy = die->uy();

  int dieArea = die->area();

  int sumTotalInstArea = 0;
  int sumMovableArea = 0;
  int sumFixedArea = 0;

  for(auto& cell : bookShelfDB_->cellVector() )
  {
    int area = cell->area();
    
    if( cell->isFixed() )
    {
      if( isOutsideDie(cell) )
        continue;
    }

    sumTotalInstArea += area;

    if( cell->isFixed() )
      sumFixedArea += area;
    else
      sumMovableArea += area;
  }

  double density = static_cast<double>(sumTotalInstArea) 
                 / static_cast<double>(dieArea);

  double util    = static_cast<double>(sumMovableArea)   
                 / static_cast<double>(dieArea - sumFixedArea);

  cout << endl;
  cout << "*** Summary of Information ***" << endl;
  cout << "---------------------------------------------" << endl;
  cout << " DESIGN INFO"                                  << endl;
  cout << "---------------------------------------------" << endl;
  cout << " DESIGN NAME      : " << name          << endl;
  cout << " NUM INSTANCE     : " << numInst       << endl;
  cout << " NUM MOVABLE      : " << numMovable()  << endl;
  cout << " NUM FIXED        : " << numFixed()    << endl;
  cout << " NUM NET          : " << numNet        << endl;
  cout << " NUM PIN          : " << numPin        << endl;
  cout << " NUM ROW          : " << numRow        << endl;
  cout << " UTIL             : " << fixed << setprecision(2) << util    * 100 << "%\n";
  cout << " DENSITIY         : " << fixed << setprecision(2) << density * 100 << "%\n";
  cout << " AREA (INSTANACE) : " << setw(10) << sumTotalInstArea << endl;
  cout << " AREA (MOVABLE)   : " << setw(10) << sumMovableArea   << endl;
  cout << " AREA (FIXED)     : " << setw(10) << sumFixedArea     << endl;
  cout << " AREA (CORE)      : " << setw(10) << dieArea          << endl;
  cout << " CORE ( " << setw(5) << dieLx  << " ";
  cout << setw(5) << dieLy  << " ) ( ";
  cout << setw(8) << dieUx  << " " << dieUy  << " )\n";
  cout << "---------------------------------------------" << endl;
}

} // namespace BookShelf
