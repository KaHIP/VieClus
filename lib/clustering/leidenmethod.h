/******************************************************************************
 * leidenmethod.h
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/


#ifndef LEIDENMETHOD_H
#define LEIDENMETHOD_H

#include "data_structure/graph_access.h"
#include "data_structure/graph_hierarchy.h"
#include "partition/partition_config.h"
#include "clustering/leiden_config.h"

#include <vector>
#include <queue>

class LeidenMethod
{
    public:
        LeidenMethod();
        virtual ~LeidenMethod();

        PartitionID performClustering(PartitionConfig &config,
                                      graph_access *G, bool = true);

        PartitionID performClusteringWithLPP(const PartitionConfig &config,
                                             graph_access *G, bool = true);

        void setLeidenConfig(const LeidenConfig &lc) { m_leiden_config = lc; }

    protected:
        void initializeSingletonClusters();

        // Phase 1: fast queue-based local moving
        NodeID performFastLocalMove(const PartitionConfig &config);

        // Phase 2: refinement within communities from phase 1
        void performRefinement(const PartitionConfig &config,
                               const std::vector<PartitionID> &flatPartition);

        graph_access *m_G;
        LeidenConfig m_leiden_config;

    private:
};

#endif // LEIDENMETHOD_H
