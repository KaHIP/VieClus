/******************************************************************************
 * evo_clustering.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#ifdef USE_MPI
#include <mpi.h>
#else
#include "tools/pseudo_mpi.h"
#endif
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "algorithms/cycle_search.h"
#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parallel_mh_clustering/parallel_mh_async_clustering.h"
#include "parse_parameters.h"
#include "clustering/leiden_config.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

int main(int argn, char **argv) {

        MPI_Init(&argn, &argv);

#ifdef USE_PSEUDO_MPI
        std::cout << "Warning: Running without MPI support (using pseudo_mpi)" << std::endl;
#endif    /* starts MPI */

        PartitionConfig partition_config;
        LeidenConfig leiden_config;
        std::string graph_filename;
        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;

        int ret_code = parse_parameters(argn, argv,
                        partition_config, graph_filename,
                        is_graph_weighted, suppress_output,
                        recursive, leiden_config); 

        if(ret_code) {
                return 0;
        }

        graph_access G;     

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);

        std::cout << "io time: " << t.elapsed()  << std::endl;
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        t.restart();
        
        partition_config.k = 1;
        parallel_mh_async_clustering mh(leiden_config);
        mh.perform_partitioning(partition_config, G);

        int rank, size;
        MPI_Comm communicator = MPI_COMM_WORLD; 
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        if( rank == ROOT ) {
                std::cout << "time spent " << t.elapsed()  << std::endl;
                G.set_partition_count(G.get_partition_count_compute());
                std::cout << "modularity \t\t\t" << ModularityMetric::computeModularity(G) << std::endl;
                
                // write the partition to the disc 
                std::stringstream filename;
                if(!partition_config.filename_output.compare("")) {
                        // no output filename given
                        filename << "tmpclustering";
                } else {
                        filename << partition_config.filename_output;
                }

                graph_io::writePartition(G, filename.str());
        }

        MPI_Finalize();
}
