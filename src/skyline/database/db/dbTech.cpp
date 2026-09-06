#include <iostream>
#include <cassert>
#include <limits>
#include <algorithm>
#include <regex>

#include "dbTech.h"

#include "dbUtil.h"
#include "dbRect.h"
#include "dbPolygon.h"

namespace db
{

dbTech::dbTech(std::shared_ptr<dbTypes> types)
  : types_               (types),
    dbu_                 (1),
    left_bus_delimiter_  ('['),
    right_bus_delimiter_ (']'),
    divider_             ('/')
{
}

dbTech::~dbTech()
{
  str2dbLayer_.clear();
  str2dbSite_.clear();
  str2dbMacro_.clear();
  str2dbViaMaster_.clear();
  str2dbViaRule_.clear();

  layers_.clear();
  sites_.clear();
  macros_.clear();
  vias_.clear();
  viaRules_.clear();
}

dbLayer*
dbTech::getLayerByName(const std::string& name)
{
  auto itr = str2dbLayer_.find(name);

  if(itr == str2dbLayer_.end())
    return nullptr;
  else
    return itr->second;
}

dbViaMaster*
dbTech::getViaMasterByName(std::string_view name)
{
  // FIXME : overhead of making std::string everytime?
  auto itr = str2dbViaMaster_.find(std::string(name));

  if(itr == str2dbViaMaster_.end())
    return nullptr;
  else
    return itr->second;
}

dbViaRule*
dbTech::getViaRuleByName(std::string_view name)
{
  // FIXME : overhead of making std::string everytime?
  auto itr = str2dbViaRule_.find(std::string(name));

  if(itr == str2dbViaRule_.end())
    return nullptr;
  else
    return itr->second;
}

dbSite*
dbTech::getSiteByName(const std::string& name)
{
  auto itr = str2dbSite_.find(name);

  if(itr == str2dbSite_.end())
    return nullptr;
  else
    return itr->second;
}

dbMacro*
dbTech::getMacroByName(const std::string& name)
{
  auto itr = str2dbMacro_.find(name);

  if(itr == str2dbMacro_.end())
    return nullptr;
  else
    return itr->second;
}

int
dbTech::getDbuLength(double micron) const
{
  double dbuLength = micron * static_cast<double>(dbu_);
  return static_cast<int>(dbuLength);
}

int
dbTech::getDbuArea(double micron) const
{
  double dbuArea = micron * static_cast<double>(dbu_) * static_cast<double>(dbu_);
  return static_cast<int>(dbuArea);
}

void
dbTech::setUnits(const lefiUnits* unit)
{
  if(unit->hasDatabase())
    dbu_ = static_cast<int>(unit->databaseNumber());
}

void
dbTech::setBusBit(const char* busBit)
{
  left_bus_delimiter_  = busBit[0];
  right_bus_delimiter_ = busBit[1];
}

void
dbTech::setDivider(const char div)
{
  divider_ = div;
}

void
dbTech::addNewLayerBookShelf(dbLayer* new_layer)
{
  int curLayerNum = layers_.size();
  new_layer->setIndex(curLayerNum);

  layers_.push_back(new_layer);

  str2dbLayer_[new_layer->name()] = new_layer;
}

void
dbTech::createNewLayer(const lefiLayer* la)
{
  if(dbu_ == 0)
  {
    std::cout << "Database Unit is not defined!" << std::endl;
    exit(1);
  }

  dbLayer* newLayer = new dbLayer;

  // Index will be used to identify vertical order of layers.
  // e.g. layer that has 0 index is the bottom layer.
  // NOTE : We assume layer starts from the bottom in the lef file.
  int curLayerNum = layers_.size();
  newLayer->setIndex(curLayerNum);

  layers_.push_back( newLayer );
  newLayer->setName(la->name());

  str2dbLayer_[newLayer->name()] = newLayer;

  if(la->hasType())      
  {
    auto type = types_->getRoutingType(std::string(la->type()));
    newLayer->setType(type);
  }

  if(la->hasDirection()) 
  {
    auto dir = types_->getLayerDirection(std::string(la->direction()));
    newLayer->setDirection(dir);
  }

  if(la->hasPitch()) 
  {
    newLayer->setXPitch(getDbuLength(la->pitch()));
    newLayer->setYPitch(getDbuLength(la->pitch()));
  }
  else if(la->hasXYPitch()) 
  { 
    newLayer->setXPitch(getDbuLength(la->pitchX()));
    newLayer->setYPitch(getDbuLength(la->pitchY()));
  }

  if(la->hasOffset()) 
  {
    newLayer->setXOffset(getDbuLength(la->offset()));
    newLayer->setYOffset(getDbuLength(la->offset()));
  }
  else if(la->hasXYOffset()) 
  {
    newLayer->setXOffset(getDbuLength(la->offsetX()));
    newLayer->setYOffset(getDbuLength(la->offsetY()));
  }

  if(la->hasWidth()) 
    newLayer->setWidth(getDbuLength(la->width()));

  if(la->hasArea())
    newLayer->setArea(getDbuArea(la->area()));

  //if(la->hasSpacing()) 
  //  newLayer->setSpacing(la->spacing());
 
  //newLayer->print();
}

void
dbTech::createNewVia(const lefiVia* via)
{
  if(dbu_ == 0)
  {
    std::cout << "Database Unit is not defined!" << std::endl;
    exit(1);
  }

  dbTechVia* newVia = new dbTechVia;
  dbViaMaster* newMaster = new dbViaMaster(newVia);

  newVia->setName(via->name());

  vias_.push_back(newMaster);
  str2dbViaMaster_[std::string(newVia->name())] = newMaster;

  if(via->hasDefault())
    newVia->setDefault();

  if(via->hasResistance())
    newVia->setResistance(via->resistance());

  const int numLayer = via->numLayers();

  // According to LEF DEF REF,
  // there must be three layers for each LEF VIA.
  assert(numLayer == 3);

  for(int l = 0; l < numLayer; l++)
  {
    dbLayer* layer 
      = getLayerByName( std::string(via->layerName(l)) );

    if(layer == nullptr)
      techNotExist("Layer", via->layerName(l));

    const int numRect = via->numRects(l);
    const int numPoly = via->numPolygons(l);
    for(int rectIdx = 0; rectIdx < numRect; rectIdx++)
    {
      int lx = getDbuLength(via->xl(l, rectIdx));
      int ly = getDbuLength(via->yl(l, rectIdx));
      int ux = getDbuLength(via->xh(l, rectIdx));
      int uy = getDbuLength(via->yh(l, rectIdx));

      dbRect* newRect = new dbRect;

      newRect->setLx(lx);
      newRect->setLy(ly);
      newRect->setUx(ux);
      newRect->setUy(uy);
      newRect->setLayer(layer);

      newVia->addRect(newRect);
    }

    for(int polyIdx = 0; polyIdx < numPoly; polyIdx++)
    {
      const auto lefPolygon = via->getPolygon(l, polyIdx);
      const int numPoints = lefPolygon.numPoints;

      dbPolygon* newPoly = new dbPolygon;

      for(int k = 0; k < numPoints; k++)
      {
        int polyX = getDbuLength( lefPolygon.x[k] );
        int polyY = getDbuLength( lefPolygon.y[k] );
        newPoly->addPoint(polyX, polyY);
      }

      newPoly->setLayer(layer);

      newVia->addPolygon(newPoly);
    }

    if(layer->type() == RoutingType::CUT)
      newVia->setCutLayer(layer);
    else
    {
      if(newVia->getBotLayer() == nullptr)
        newVia->setBotLayer(layer);
      else
      {
        const auto curBotLayer = newVia->getBotLayer();
        int botIdx = curBotLayer->index();
        int topIdx = layer->index();

        if(botIdx < topIdx)
          newVia->setTopLayer(layer);
        else // botIdx > topIdx
        {
          newVia->setBotLayer(layer);
          newVia->setTopLayer(curBotLayer);
        }
      }
    }
  }

  // To find shape in the cut layer of this via.
  newVia->setBoundary();

  // Sanity Check
  assert(newVia->getBotLayer()->index() 
       < newVia->getTopLayer()->index());

  if(via->hasViaRule())
    assert(0); // Not Supported
}

void
dbTech::createNewViaRule(const lefiViaRule* viaRule)
{
  const char* via_rule_name = viaRule->name();

  dbViaRule* new_via_rule = new dbViaRule(via_rule_name);

  viaRules_.push_back(new_via_rule);
  str2dbViaRule_[std::string(new_via_rule->getName())] = new_via_rule;
}

void
dbTech::createNewSite(const lefiSite* site)
{
  if(dbu_ == 0)
  {
    std::cout << "Database Unit is not defined!" << std::endl;
    exit(1);
  }

  dbSite* newSite = new dbSite;
  sites_.push_back(newSite);

  newSite->setName(site->name());
  str2dbSite_[newSite->name()] = newSite;

  if(site->hasSize())
  {
    newSite->setSizeX( getDbuLength(site->sizeX()) );
    newSite->setSizeY( getDbuLength(site->sizeY()) );
  }

  if(site->hasClass())
  {
    auto siteClass = types_->getSiteClass(std::string(site->siteClass()));
    newSite->setSiteClass( siteClass );
  }

  if(site->hasXSymmetry())
    newSite->setSymmetryX(true);
  if(site->hasYSymmetry())
    newSite->setSymmetryY(true);
  if(site->has90Symmetry())
    newSite->setSymmetryR90(true);

  //newSite->print();
}

void
dbTech::addPinToMacro(const lefiPin* pi, dbMacro* topMacro)
{
  if(dbu_ == 0)
  {
    std::cout << "Database Unit is not defined!" << std::endl;
    exit(1);
  }

  dbMTerm* newMTerm = new dbMTerm;
  newMTerm->setName( std::string(pi->name()) );
  newMTerm->setMacro( topMacro );

  if(pi->hasDirection() )
  {
    auto dir = types_->getPinDirection( std::string(pi->direction()) );
    newMTerm->setPinDirection( dir );
  }

  if(pi->hasUse())
  {
    auto use = types_->getPinUsage( std::string(pi->use()) );
    newMTerm->setPinUsage(use);
  }

  if(pi->hasShape())
  {
    auto shape = types_->getPinShape( std::string(pi->shape()) );
    newMTerm->setPinShape(shape);
  }

  int numPorts = pi->numPorts();

  dbLayer* curLayer = nullptr;
  for(int i = 0; i < numPorts; i++)
  {
    lefiGeometries* geo = pi->port(i);
    lefiGeomRect*    lrect = nullptr;
    lefiGeomPolygon* lpoly = nullptr;

    int numItems = geo->numItems();
    for(int j = 0; j < numItems; j++)
    {
      switch(geo->itemType(j))
      {
        case lefiGeomLayerE:
        {
          curLayer = getLayerByName( std::string(geo->getLayer(j)) );

          if(curLayer == nullptr)
            techNotExist("Layer", geo->getLayer(j));
          break;
        }
        case lefiGeomRectE:
        {
          lrect = geo->getRect(j);

          int rectLx = getDbuLength( lrect->xl );
          int rectLy = getDbuLength( lrect->yl );
          int rectUx = getDbuLength( lrect->xh );
          int rectUy = getDbuLength( lrect->yh );
 
          dbMTermPort* newPort = new dbMTermPort;
          newPort->addPoint(rectLx, rectLy);
          newPort->addPoint(rectUx, rectLy);
          newPort->addPoint(rectUx, rectUy);
          newPort->addPoint(rectLx, rectUy);

          newPort->setLayer(curLayer);
          // POLYGON implicitly finish with the starting point.
          newMTerm->addPort( newPort );
          break;
        }

        case lefiGeomPolygonE:
        {
          lpoly = geo->getPolygon(j);

          const int numPoints = lpoly->numPoints;
          dbMTermPort* newPort = new dbMTermPort;

          for(int k = 0; k < numPoints; k++)
          {
            int polyX = getDbuLength( lpoly->x[k] );
            int polyY = getDbuLength( lpoly->y[k] );
            newPort->addPoint(polyX, polyY);
          }

          newMTerm->addPort( newPort );
          break;
        }

        case lefiGeomPolygonIterE:
        {
          assert(0);
          break;
        }
  
        default:
          break;
      }
    }
  }
  newMTerm->setBoundary();
  topMacro->addMTerm( newMTerm );
}

void
dbTech::addObsToMacro(const lefiObstruction* obs, dbMacro* topMacro)
{
  lefiGeometries* geo = obs->geometries();

  int numItems = geo->numItems();
  dbLayer* curLayer = nullptr;
  for(int i = 0; i < numItems; i++)
  {
    lefiGeomRect* lrect = nullptr;
    lefiGeomPolygon* lpoly = nullptr;

    switch(geo->itemType(i)) 
    {
      case lefiGeomLayerE:
      {
        curLayer = getLayerByName( std::string(geo->getLayer(i)) );

        if(curLayer == nullptr)
          techNotExist("Layer", geo->getLayer(i));
        break;
      }
      case lefiGeomRectE:
      {
        lrect = geo->getRect(i);

        int rectLx = getDbuLength( lrect->xl );
        int rectLy = getDbuLength( lrect->yl );
        int rectUx = getDbuLength( lrect->xh );
        int rectUy = getDbuLength( lrect->yh );
 
        dbObs* newObs = new dbObs; 
        newObs->addPoint(rectLx, rectLy);
        newObs->addPoint(rectUx, rectLy);
        newObs->addPoint(rectUx, rectUy);
        newObs->addPoint(rectLx, rectUy);

        newObs->setLayer(curLayer);
        // POLYGON implicitly finish with the starting point.
        topMacro->addObs( newObs );
        break;
      }

      case lefiGeomPolygonE:
      {
        lpoly = geo->getPolygon(i);

        int numPoints = lpoly->numPoints;
        dbObs* newObs = new dbObs;

        for(int k = 0; k < numPoints; k++)
        {
          int polyX = getDbuLength( lpoly->x[k] );
          int polyY = getDbuLength( lpoly->y[k] );
          newObs->addPoint(polyX, polyY);
        }

        topMacro->addObs( newObs );
        break;
      }

      case lefiGeomPolygonIterE:
      {
        assert(0);
        break;
      }

      default: 
        break;    
    }
  }
}

void 
dbTech::addEdgeSpacing(int type1, int type2, double spacing)
{
  edge_spacing_.push_back({type1, type2, spacing});
}

dbMacro*
dbTech::getNewMacro(const char* name)
{
  if(dbu_ == 0)
  {
    std::cout << "Database Unit is not defined!" << std::endl;
    exit(1);
  }

  dbMacro* new_macro = new dbMacro;
  macros_.push_back(new_macro);

  new_macro->setName(name);
  str2dbMacro_[new_macro->name()] = new_macro;

  return new_macro;
}

void
dbTech::fillNewMacro(const lefiMacro* ma, dbMacro* new_macro)
{
  if( ma->hasClass() )
  {
    auto macroClass = types_->getMacroClass(std::string(ma->macroClass()));
    new_macro->setMacroClass(macroClass);
  }

  if( ma->hasOrigin() ) 
  {
    new_macro->setOrigX( getDbuLength(ma->originX()) );
    new_macro->setOrigY( getDbuLength(ma->originY()) );
  }

  if( ma->hasSize() )
  {
    new_macro->setSizeX( getDbuLength(ma->sizeX()) );
    new_macro->setSizeY( getDbuLength(ma->sizeY()) );
  }

  if( ma->hasSiteName() ) 
  {
    dbSite* site = getSiteByName( std::string(ma->siteName()) );
    new_macro->setSite( site );
  }

  if( ma->hasXSymmetry()  ) new_macro->setSymmetryX(true);
  if( ma->hasYSymmetry()  ) new_macro->setSymmetryY(true);
  if( ma->has90Symmetry() ) new_macro->setSymmetryR90(true);

  const int num_properties = ma->numProperties();
  if( num_properties > 0 ) 
  {
    for(int i = 0; i < num_properties; i++) 
    {
      std::string prop_name(ma->propName(i));
      if(prop_name == "LEF58_EDGETYPE")
      {
        std::string s(ma->propValue(i));
        std::regex rgx_left("LEFT (\\d+) ;");
        std::regex rgx_right("RIGHT (\\d+) ;");
  
        std::smatch match;
        int left_edge_type, right_edge_type;
        if(std::regex_search(s, match, rgx_left) && match.size() > 1) 
          std::istringstream(match.str(1)) >> left_edge_type;
        if(std::regex_search(s, match, rgx_right) && match.size() > 1) 
          std::istringstream(match.str(1)) >> right_edge_type;
        
        new_macro->setLEdgeType(left_edge_type);
        new_macro->setREdgeType(right_edge_type);
      }
    }
  }
  //new_macro->print();
}
  
}
