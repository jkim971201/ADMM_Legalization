#ifndef SKYLINE_CHRONO_H
#define SKYLINE_CHRONO_H

#include <chrono>

namespace util
{

std::chrono::high_resolution_clock::time_point getChronoNow();

double evalTime(const std::chrono::system_clock::time_point& start);

double evalTimeMilliSeconds(const std::chrono::system_clock::time_point& start);

}

#endif
