#pragma once
#include "utils.cuh"
#include "modularity_optimisation.cuh"

namespace Louvain
{

int detectSparseCluster(
  int V,
  int minNumVertex, 
  HostStructure& hostStructures, 
  DeviceStructure& deviceStructures);

bool removeSparseCluster(
  float minGain, 
	int minSize,
  DeviceStructure& deviceStructures, 
  HostStructure& hostStructures);
};
