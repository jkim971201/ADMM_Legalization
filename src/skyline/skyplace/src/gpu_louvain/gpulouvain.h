#pragma once

#include "utils.cuh"
#include "modularity_optimisation.cuh"
#include "community_aggregation.cuh"

namespace Louvain
{

struct HostStructure;

void gpuLouvain(
  HostStructure& hostStructures, 
	float minGain, 
	int minClusterSize = 3);

};
