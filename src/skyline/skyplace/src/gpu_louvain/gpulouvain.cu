#include <stdio.h>
#include "gpulouvain.h"
#include "utils.cuh"
#include "modularity_optimisation.cuh"
#include "community_aggregation.cuh"
#include "removeSparse.cuh"

namespace Louvain
{

void gpuLouvain(HostStructure& hostStructures, float minGain, int minClusterSize)
{
  DeviceStructure deviceStructures;
  AggregationStructure aggregationPhaseStructures;

  copyStructures(hostStructures, deviceStructures, aggregationPhaseStructures);
  initM(hostStructures);

  int i = 0;
  int maxIter = 500;
  for(;;) 
  {
    bool worthContinue = optimiseModularity(minGain, deviceStructures, hostStructures);

		if(!worthContinue)
			break;

    aggregateCommunities(deviceStructures, hostStructures, aggregationPhaseStructures);

    int V;
    HANDLE_ERROR(cudaMemcpy(&V, deviceStructures.V, sizeof(int), cudaMemcpyDeviceToHost));
    //printf("Louvain Iter[%02d] V: %d Modularity: %f\n", i, V, calculateModularity(V, hostStructures.M, deviceStructures));
    i++;
    if(i == maxIter)
      break;
  }

  int V;
  HANDLE_ERROR(cudaMemcpy(&V, deviceStructures.V, sizeof(int), cudaMemcpyDeviceToHost));

  copyDS2HS(deviceStructures, hostStructures);
  deleteStructures(deviceStructures, aggregationPhaseStructures);
}

};
