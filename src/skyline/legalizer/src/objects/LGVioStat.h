#ifndef LG_VIO_STAT_H
#define LG_VIO_STAT_H

namespace legalizer
{

struct LGVioStat
{
  int total_ovf;
  int num_row_orient_vio;
  int num_pin_short_vio;
  int num_pin_access_vio;
  int num_edge_spacing_vio;
};

}

#endif
