#include "util/Chrono.h"

namespace util
{

std::chrono::high_resolution_clock::time_point getChronoNow()
{
  return std::chrono::high_resolution_clock::now();
}

double evalTime(const std::chrono::system_clock::time_point& start)
{
  auto time_here = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration_double = time_here - start;
  return duration_double.count();
}

double evalTimeMilliSeconds(const std::chrono::system_clock::time_point& start)
{
  auto time_here = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration_double = time_here - start;
  return duration_double.count();
}

}
