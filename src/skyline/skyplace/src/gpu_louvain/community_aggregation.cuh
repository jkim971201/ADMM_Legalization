#ifndef LOUVIAN_COMMUNITY_AGGREGATION_CUH
#define LOUVIAN_COMMUNITY_AGGREGATION_CUH

#include "utils.cuh"

namespace Louvain
{

/**
 * Transforms every community into single vertex.
 * @param deviceStructures structures kept in device memory
 * @param hostStructures   structures kept in host memory
 */
void aggregateCommunities(
  DeviceStructure&      deviceStructures, 
  HostStructure&        hostStructures,
	AggregationStructure& aggregationPhaseStructures);
};

#endif //LOUVIAN_COMMUNITY_AGGREGATION_CUH
