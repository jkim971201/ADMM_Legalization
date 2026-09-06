#ifndef SKYPLACE_UTIL_H
#define SKYPLACE_UTIL_H

#include <chrono>

namespace skyplace
{

std::chrono::high_resolution_clock::time_point getChronoNow();

double evalTime(const std::chrono::system_clock::time_point& start);


}

#endif
