#include <cassert>
#include <cstdio>

#include "skyline/SkyLine.h"

#include "db/dbDatabase.h"
#include "skygui/SkyGui.h"
#include "skyplace/SkyPlace.h"
#include "legalizer/Legalizer.h"

namespace skyline
{

// Singleton Design Pattern
SkyLine* SkyLine::getStaticPtr()
{
  static SkyLine sky;
  return &sky;
}

SkyLine::SkyLine()
{
  db_        = std::make_shared<db::dbDatabase>();
  gui_       = std::make_unique<gui::SkyGui>(db_);
  skyplace_  = std::make_unique<skyplace::SkyPlace>(db_);
  legalizer_ = std::make_unique<legalizer::Legalizer>(db_);
}

SkyLine::~SkyLine() {}

void
SkyLine::readLef(const char* file_path)
{
  db_->readLef(file_path);
}

void
SkyLine::readDef(const char* file_path)
{
  db_->readDef(file_path);
}

void
SkyLine::readVerilog(const char* file_path)
{
  db_->readVerilog(file_path);
}

void
SkyLine::readBookShelf(const char* file_path)
{
  db_->readBookShelf(file_path);
}

void
SkyLine::writeDef(const char* file_path)
{
  db_->writeDef(file_path);
}

void
SkyLine::writeBookShelf(const char* file_path)
{
  db_->writeBookShelf(file_path);
}

void
SkyLine::setTopModuleName(const char* top_name)
{
  db_->setTopModuleName(top_name);
}

skyplace::SkyPlace*
SkyLine::getSkyPlaceGP()
{
  return skyplace_.get();
}

legalizer::Legalizer*
SkyLine::getSkyPlaceLG()
{
  return legalizer_.get();
}

void
SkyLine::display()
{
  gui_->display();
}

} // namespace skyline
