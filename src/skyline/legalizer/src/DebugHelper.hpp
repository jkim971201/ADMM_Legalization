#include "cuda_linalg/CudaVectorAlgebra.h"

namespace legalizer
{

void printPixelUsage(
  const int x_size,
  const int y_size,
  const cuda_linalg::CudaVector<int>& pixel_usage)
{
  std::vector<std::pair<int, int>> minus;

  std::vector<int> host_usage(pixel_usage.size());

  thrust::copy(pixel_usage.begin(), pixel_usage.end(), host_usage.begin());

  printf("Pixel Usage %d X %d\n", x_size, y_size);
  for(int y = y_size - 1; y >= 0; y--)
  {
    for(int x = 0; x < x_size; x++)
    {
      int u = host_usage[x * y_size + y];
      if(u < 0)
      {
        printf("X");
        minus.push_back({x, y});
      }
      else
        printf("%d", u);
    }
    printf("\n");
  }

  if(minus.empty() == false)
  {
    printf("Minus found!\n");
    for(const auto& [x, y] : minus)
      printf("(%d, %d)\n", x, y);
  }
}

}
