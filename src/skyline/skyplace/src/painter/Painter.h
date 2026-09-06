#ifndef PAINTER_H
#define PAINTER_H

#include <filesystem>
#include <unordered_map>
#include <QColor>

class QPainter;
class QApplication;

namespace skyplace 
{

class SkyPlaceDB;

class Painter 
{
  public:

    Painter(std::shared_ptr<SkyPlaceDB> db);
    ~Painter();

    void preparePlot(const std::string& path);
    void saveImage(int iter, float hpwl, float overflow);

  private:

    void drawDieRect(float scale, QPainter* painter) const;

    void drawFixedRect(float k_scale, QPainter* painter) const;

    void drawInstRect(float k_scale, QPainter* painter) const;

    std::shared_ptr<SkyPlaceDB> db_;

    std::unique_ptr<QApplication> qapp_;

    std::filesystem::path path_to_save_image_;

    std::unordered_map<int, QColor> clusterId2Color_;
};

} // namespace skyline

#endif
