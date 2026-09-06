#ifndef SKYLINE_H
#define SKYLINE_H

#include <memory>
#include <string_view>

namespace db { 
  class dbDatabase; 
}

namespace gui { 
  class SkyGui;    
}

namespace skyplace { 
  class SkyPlace;   
}

namespace legalizer
{
  class Legalizer;
}

namespace skyline
{

class SkyLine
{
  public:

    static SkyLine* getStaticPtr();

    void readLef          (const char* file_path);
    void readDef          (const char* file_path);
    void readVerilog      (const char* file_path);
    void readBookShelf    (const char* file_path);

    void writeDef         (const char* file_path);
    void writeBookShelf   (const char* file_path);

    void setTopModuleName (const char*  top_name);

    // Global Placer
    skyplace::SkyPlace* getSkyPlaceGP();

    // Legalizer
    legalizer::Legalizer* getSkyPlaceLG();

    void display();

  private:

    // For Singleton, we should make Constructor and Deconstructor as private
    SkyLine();
    ~SkyLine();

    std::shared_ptr<db::dbDatabase>       db_;
    std::unique_ptr<gui::SkyGui>          gui_;
    std::unique_ptr<legalizer::Legalizer> legalizer_;
    std::unique_ptr<skyplace::SkyPlace>   skyplace_;
};

}

#endif
