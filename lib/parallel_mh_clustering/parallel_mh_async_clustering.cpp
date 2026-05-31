/******************************************************************************
 * parallel_mh_async_clustering.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/

#include <algorithm>
#include <fstream>
#include <iostream>
#ifdef USE_MPI
#include <mpi.h>
#else
#include "tools/pseudo_mpi.h"
#endif
#include <sstream>
#include <stdio.h>

#include "diversifyer.h"
#include "exchange/exchanger_clustering.h"
#include "graph_io.h"
#include "graph_partitioner.h"
#include "parallel_mh_async_clustering.h"
#include "quality_metrics.h"
#include "random_functions.h"

parallel_mh_async_clustering::parallel_mh_async_clustering() : MASTER(0), m_time_limit(0) {
        m_best_global_objective = std::numeric_limits<EdgeWeight>::max();
        m_best_cycle_objective  = std::numeric_limits<EdgeWeight>::max();
        m_rounds                = 0;
        m_termination           = false;
        m_communicator          = MPI_COMM_WORLD;
        MPI_Comm_rank( m_communicator, &m_rank);
        MPI_Comm_size( m_communicator, &m_size);
}

parallel_mh_async_clustering::parallel_mh_async_clustering(const LeidenConfig & leiden_config) : MASTER(0), m_time_limit(0), m_leiden_config(leiden_config) {
        m_best_global_objective = std::numeric_limits<EdgeWeight>::max();
        m_best_cycle_objective  = std::numeric_limits<EdgeWeight>::max();
        m_rounds                = 0;
        m_termination           = false;
        m_communicator          = MPI_COMM_WORLD;
        MPI_Comm_rank( m_communicator, &m_rank);
        MPI_Comm_size( m_communicator, &m_size);
}

parallel_mh_async_clustering::parallel_mh_async_clustering(MPI_Comm communicator) : MASTER(0), m_time_limit(0) {
        m_best_global_objective = std::numeric_limits<EdgeWeight>::max();
        m_best_cycle_objective  = std::numeric_limits<EdgeWeight>::max();
        m_rounds                = 0;
        m_termination           = false;
        m_communicator          = communicator;
        MPI_Comm_rank( m_communicator, &m_rank);
        MPI_Comm_size( m_communicator, &m_size);

}

parallel_mh_async_clustering::~parallel_mh_async_clustering() {
        delete[] m_best_global_map;
}

void parallel_mh_async_clustering::perform_partitioning(const PartitionConfig & partition_config, graph_access & G) {
        m_time_limit      = partition_config.time_limit;
        m_island          = new population_clustering(m_communicator, partition_config, m_leiden_config);
        m_best_global_map = new PartitionID[G.number_of_nodes()];

        srand(partition_config.seed*m_size+m_rank);
        random_functions::setSeed(partition_config.seed*m_size+m_rank);

        PartitionConfig ini_working_config  = partition_config; 
        initialize( ini_working_config, G);

        global_timer_restart();
        exchanger_clustering ex(m_communicator);
        do {
                PartitionConfig working_config  = partition_config; 

                perform_local_partitioning( working_config, G );

                //push and recv 
                if( global_timer_elapsed() <= m_time_limit && m_size > 1) {
                        unsigned messages = ceil(log(m_size));
                        for( unsigned i = 0; i < messages; i++) {
                                ex.push_best( working_config, G, *m_island );
                                ex.recv_incoming( working_config, G, *m_island );
                        }
                }

                m_rounds++;
        } while( global_timer_elapsed() <= m_time_limit );

        collect_best_partitioning(G, partition_config);


        m_island->print();

        //print logfile (for convergence plots)
        if( partition_config.mh_print_log ) {
                std::stringstream filename_stream;
                filename_stream << "log_"<<  partition_config.graph_filename <<   
                        "_m_rank_" <<  m_rank <<  
                        "_file_" <<  
                        "_seed_" <<  partition_config.seed <<  
                        "_k_" <<  partition_config.k;

                std::string filename(filename_stream.str());
                m_island->write_log(filename);
        }

        delete m_island;
}

void parallel_mh_async_clustering::initialize(PartitionConfig & working_config, graph_access & G) {
        // each PE performs a clustering 
        // estimate the runtime of a partitioner call 
        // calculate the poolsize and async Bcast the poolsize.
        // recv. has to be sync
        Individuum first_one;
        global_timer_restart();
        m_island->createIndividuum( working_config, G, first_one, true); 
        std::cout <<  "created with objective " <<  first_one.objective << std::endl;

        double time_spend = global_timer_elapsed();
        m_island->insert(G, first_one);

        //compute S and Bcast
        int population_size = 1;
        double fraction     = working_config.mh_initial_population_fraction;
        int POPSIZE_TAG     = 10;

        if( m_rank == ROOT ) {
                double fraction_to_spend_for_IP = (double)m_time_limit / fraction;
                population_size                 = ceil(fraction_to_spend_for_IP / time_spend);

                for( int target = 1; target < m_size; target++) {
                        MPI_Request rq;
                        MPI_Isend(&population_size, 1, MPI_INT, target, POPSIZE_TAG, m_communicator, &rq); 
                }
        } else {
                MPI_Status rst;
                MPI_Recv(&population_size, 1, MPI_INT, ROOT, POPSIZE_TAG, m_communicator, &rst); 
        }

        population_size = std::min(std::max(10, population_size),(int)working_config.mh_pool_size);
        std::cout <<  "poolsize = " <<  population_size  << std::endl;

        //set S
        m_island->set_pool_size(population_size);
        working_config.mh_pool_size = population_size;

}

double parallel_mh_async_clustering::collect_best_partitioning(graph_access & G, const PartitionConfig & config) {
        //perform partitioning locally
        double max_objective = 0;
        m_island->apply_fittest(G, max_objective);

        double best_local_objective   = max_objective;
        double best_local_objective_m = max_objective;
        double best_global_objective  = 0;

        PartitionID* best_local_map = new PartitionID[G.number_of_nodes()];

        forall_nodes(G, node) {
                best_local_map[node] = G.getPartitionIndex(node);
        } endfor

        MPI_Allreduce(&best_local_objective_m, &best_global_objective, 1, MPI_DOUBLE, MPI_MAX, m_communicator);

        // now we know what the best objective is ... find the best balance
        int bcaster = best_local_objective == best_global_objective ? m_rank : std::numeric_limits<int>::max();
        int g_bcaster = 0;

        MPI_Allreduce(&bcaster, &g_bcaster, 1, MPI_INT, MPI_MIN, m_communicator);
        MPI_Bcast(best_local_map, G.number_of_nodes(), MPI_INT, g_bcaster, m_communicator);

        forall_nodes(G, node) {
                G.setPartitionIndex(node, best_local_map[node]);
        } endfor

        delete[] best_local_map;

        return best_global_objective;
}

double parallel_mh_async_clustering::perform_local_partitioning(PartitionConfig & working_config, graph_access & G) {

        quality_metrics qm;
        unsigned local_repetitions = working_config.local_partitioning_repetitions;

        //start a new round
        for( unsigned i = 0; i < local_repetitions; i++) {
                if( !m_island->is_full()) {
                        Individuum first_ind;
                        m_island->createIndividuum(working_config, G, first_ind, true);
                        m_island->insert(G, first_ind);
                } else {
                        //perform combine operations
                        Individuum first_rnd;
                        Individuum second_rnd;
                        Individuum output;
                        m_island->get_two_individuals_tournament(first_rnd, second_rnd);
                        
                        int decision = random_functions::nextInt(0,86);
                        if( 0 <= decision && decision <= 20) {
                              m_island->combine_basic_flat(working_config, G, first_rnd, second_rnd, output);
                        }        
                        if( 21 <= decision && decision <= 40) {
                              m_island->combine_improved_flat(working_config, G, first_rnd, second_rnd, output);
                        }
                        if( 41 <= decision && decision <= 60) {
                              m_island->combine_improved_flat_with_sclp(working_config, G, first_rnd, output);
                        } 
                        if( 61 <= decision && decision <= 80) {
                              m_island->combine_improved_multilevel(working_config, G, first_rnd, second_rnd, output);
                        }
                        if( 81 <= decision && decision <= 83) {
                              m_island->combine_improved_flat_with_partitioning(working_config, G, first_rnd, output);
                        } 
                        if( 84 <= decision && decision <= 86) {
                              m_island->mutate(working_config, G, first_rnd, second_rnd, output);
                        }

                        
                        m_island->insert(G, output);
                }

                //try to combine to random inidividuals from pool 
                if( global_timer_elapsed() > m_time_limit ) {
                        break;
                }
        }

        double max_objective = 0;
        m_island->apply_fittest(G, max_objective);

        return max_objective;
}


