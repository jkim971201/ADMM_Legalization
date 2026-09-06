#ifndef DB_RECT_H
#define DB_RECT_H

#include "dbBox.h"

namespace db
{

class dbLayer;

// This class is used when
// #1. for describing the pin shape of dbMTerm (LEF PIN)
// #2. for describing the obstruct shape       (LEF OBS)
// #3. for describing the boundary of dbBTerm  (DEF PIN PORT)

// dbRect has "dbLayer", while dbBox does not.
class dbRect : public dbBox
{
  public:

  dbRect()
    : dbBox(0, 0, 0, 0), layer_ (nullptr)
  {}

  dbRect(int lx, int ly, int ux, int uy, dbLayer* _layer)
    : dbBox(lx, ly, ux, uy), layer_ (_layer) 
  {}

  // Setters
  void setLayer(dbLayer* layer) { layer_ = layer; }

  // Getters
  const dbLayer* layer() const { return layer_; }
        dbLayer* layer()       { return layer_; }

  protected:

    dbLayer* layer_;
};

}

#endif
