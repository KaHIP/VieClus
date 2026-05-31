/******************************************************************************
 * vieclus_interface.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/

#include <iostream>
#include <limits>

#include "vieclus_interface.h"
#include "data_structure/graph_access.h"
#include "tools/random_functions.h"
#include "tools/modularitymetric.h"
#include "tools/global_timer.h"
#include "parallel_mh_clustering/parallel_mh_async_clustering.h"
#include "clustering/leiden_config.h"
#include "partition/partition_config.h"
#include "macros_assertions.h"
#include "configuration.h"

#ifdef USE_MPI
#include <mpi.h>
#else
#include "tools/pseudo_mpi.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

static void internal_build_graph(int* n,
                                 int* vwgt,
                                 int* xadj,
                                 int* adjcwgt,
                                 int* adjncy,
                                 graph_access & G) {
        G.build_from_metis(*n, xadj, adjncy);

        if (vwgt != NULL) {
                forall_nodes(G, node) {
                        G.setNodeWeight(node, vwgt[node]);
                } endfor
        }

        if (adjcwgt != NULL) {
                forall_edges(G, e) {
                        G.setEdgeWeight(e, adjcwgt[e]);
                } endfor
        }
}

void vieclus_clustering(int* n, int* vwgt, int* xadj,
                        int* adjcwgt, int* adjncy,
                        bool suppress_output, int seed,
                        double time_limit, int cluster_upperbound,
                        bool leiden, double leiden_theta,
                        double* modularity, int* num_clusters, int* clustering) {

        int argn_dummy = 0;
        char** argv_dummy = NULL;
        MPI_Init(&argn_dummy, &argv_dummy);

#ifdef _OPENMP
        omp_set_num_threads(1);
#endif

        // Configure (always use STRONG)
        PartitionConfig config;
        configuration cfg;
        cfg.strong(config);

        config.k = 1;
        config.seed = seed;
        config.time_limit = time_limit;
        config.suppress_partitioner_output = suppress_output;

        if (cluster_upperbound > 0) {
                config.cluster_upperbound = cluster_upperbound;
                config.upper_bound_partition = cluster_upperbound;
        } else {
                config.cluster_upperbound = std::numeric_limits<NodeWeight>::max() / 2;
                config.upper_bound_partition = std::numeric_limits<NodeWeight>::max() / 2;
        }

        // Build graph
        graph_access G;
        internal_build_graph(n, vwgt, xadj, adjcwgt, adjncy, G);
        G.set_partition_count(1);

        srand(seed);
        random_functions::setSeed(seed);

        // Run clustering
        global_timer_restart();
        LeidenConfig leiden_config(leiden, leiden_theta);
        parallel_mh_async_clustering mh(leiden_config);
        mh.perform_partitioning(config, G);

        // Compute results
        G.set_partition_count(G.get_partition_count_compute());
        *modularity = ModularityMetric::computeModularity(G);
        *num_clusters = G.get_partition_count();

        forall_nodes(G, node) {
                clustering[node] = G.getPartitionIndex(node);
        } endfor

        MPI_Finalize();
}
