#include "dbBTerm.h"

#include <limits>
#include <iostream>
#include <cassert>

namespace db
{

void
dbBTermPort::setLocation()
{
  switch(orient_)
  {
    case Orient::N :
    {
      lx_ = origX_ + offsetLx_;
      ly_ = origY_ + offsetLy_;
      ux_ = origX_ + offsetUx_;
      uy_ = origY_ + offsetUy_;
      break;
    }
    case Orient::S :
    {
      lx_ = origX_ - offsetUx_;
      ly_ = origY_ - offsetUy_;
      ux_ = origX_ - offsetLx_;
      uy_ = origY_ - offsetLy_;
      break;
    }
    case Orient::W :
    {
      lx_ = origX_ - offsetUy_;
      ly_ = origY_ + offsetLx_;
      ux_ = origX_ - offsetLy_;
      uy_ = origY_ + offsetUx_;
      break;
    }
    case Orient::E :
    {
      lx_ = origX_ + offsetLy_;
      ly_ = origY_ - offsetUx_;
      ux_ = origX_ + offsetUy_;
      uy_ = origY_ - offsetLx_;
      break;
    }
    case Orient::FN :
    {
      assert(0);
      break;
    }
    case Orient::FS :
    {
      assert(0);
      break;
    }
    case Orient::FW :
    {
      assert(0);
      break;
    }
    case Orient::FE :
    {
      assert(0);
      break;
    }
  } 
}

dbBTerm::dbBTerm()
 : lx_       (std::numeric_limits<int>::max()),
   ly_       (std::numeric_limits<int>::max()),
   ux_       (std::numeric_limits<int>::min()),
   uy_       (std::numeric_limits<int>::min()),
   name_     (""),
   net_      (nullptr),
   direction_(PinDirection::INOUT)
{}

void
dbBTerm::addPort(dbBTermPort* port) 
{ 
  lx_ = std::min(lx_, port->lx());
  ly_ = std::min(ly_, port->ly());
  ux_ = std::max(ux_, port->ux());
  uy_ = std::max(uy_, port->uy());

  port->setBTerm(this); 
  ports_.push_back(port); 
}

void
dbBTerm::print() const
{
  std::cout << std::endl;
  std::cout << "PIN NAME  : " << name_ << std::endl;
  std::cout << "NET NAME  : " << net_->name() << std::endl;
  std::cout << "LAYER     : " << ports_[0]->layer()->name() << std::endl;
  std::cout << "ORIG X    : " << ports_[0]->origX() << std::endl;
  std::cout << "ORIG Y    : " << ports_[0]->origY() << std::endl;
  std::cout << "OFFSET LX : " << ports_[0]->offsetLx() << std::endl;
  std::cout << "OFFSET LY : " << ports_[0]->offsetLy() << std::endl;
  std::cout << "OFFSET UX : " << ports_[0]->offsetUx() << std::endl;
  std::cout << "OFFSET UY : " << ports_[0]->offsetUy() << std::endl;
  std::cout << "LX        : " << lx() << std::endl;
  std::cout << "LY        : " << ly() << std::endl;
  std::cout << "UX        : " << ux() << std::endl;
  std::cout << "UY        : " << uy() << std::endl;
  std::cout << "CX        : " << cx() << std::endl;
  std::cout << "CY        : " << cy() << std::endl;
  std::cout << std::endl;
}

}
