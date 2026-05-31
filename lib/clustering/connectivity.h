/******************************************************************************
 * connectivity.h
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/


#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#include "data_structure/graph_access.h"

// Splits disconnected communities into connected components.
// Modifies partition indices in-place. Returns new partition count.
PartitionID splitDisconnectedCommunities(graph_access &G);

#endif // CONNECTIVITY_H
