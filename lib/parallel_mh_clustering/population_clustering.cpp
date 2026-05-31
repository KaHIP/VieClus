/******************************************************************************
 * population_clustering.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#ifdef USE_MPI
#include <mpi.h>
#else
#include "tools/pseudo_mpi.h"
#endif
#include <sstream>
#include <set>

#include "diversifyer.h"
#include "balance_configuration.h"
#include "graph_partitioner.h"
#include "population_clustering.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "tools/global_timer.h"
#include "uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "clustering/louvainmethod.h"
#include "clustering/leidenmethod.h"
#include "clustering/connectivity.h"
#include "clustering/coarsening/contractor.h"
#include "partition/coarsening/clustering/size_constraint_label_propagation.h"
#include "tools/modularitymetric.h"
#include "tools/graph_extractor.h"

population_clustering::population_clustering( MPI_Comm communicator, const PartitionConfig & partition_config, const LeidenConfig & leiden_config ) {
        m_population_clustering_size    = partition_config.mh_pool_size;
        m_no_partition_calls = 0;
        m_num_NCs            = partition_config.mh_num_ncs_to_compute;
        m_num_NCs_computed   = 0;
        m_num_ENCs           = 0;
        m_time_stamp         = 0;
        m_communicator       = communicator;
        m_use_leiden         = leiden_config.enabled;
        m_leiden_theta       = leiden_config.theta;
        global_timer_restart();
        best_objective = -1;
}

void population_clustering::do_leiden_clustering(PartitionConfig & config, graph_access & G, bool start_w_singletons) {
        LeidenMethod lm;
        lm.setLeidenConfig(LeidenConfig{true, m_leiden_theta});
        lm.performClustering(config, &G, start_w_singletons);
}

population_clustering::~population_clustering() {
        for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                delete[] (m_internal_population_clustering[i].partition_map);
                delete m_internal_population_clustering[i].cut_edges;
        }         
}

void population_clustering::set_pool_size(int size) {
        m_population_clustering_size = size;
}

void population_clustering::createIndividuum(const PartitionConfig & config, 
                graph_access & G, 
                Individuum & ind, bool output) {

        PartitionConfig copy = config;

        copy.upper_bound_partition                     = G.number_of_nodes() + 1;
        copy.graph_allready_partitioned                = false;
        copy.node_ordering                             = RANDOM_NODEORDERING;
        copy.lm_cluster_coarsening_factor              = 1;

        copy.cluster_upperbound = random_functions::nextInt(G.number_of_nodes()/10, G.number_of_nodes());
        copy.upper_bound_partition = copy.cluster_upperbound; 
        copy.cluster_coarsening_factor = 1;

        int lp_levels = random_functions::nextInt(0,10);
        if( 0 <= lp_levels && lp_levels <= 7) 
                copy.lm_number_of_label_propagation_levels = 0; 
        if( lp_levels == 8) 
                copy.lm_number_of_label_propagation_levels = 1; 
        if( lp_levels == 9) 
                copy.lm_number_of_label_propagation_levels = 2; 
        if( lp_levels == 10) 
                copy.lm_number_of_label_propagation_levels = 3; 

        if (m_use_leiden) {
                LeidenMethod lm;
                lm.setLeidenConfig(LeidenConfig{true, m_leiden_theta});
                lm.performClustering(copy, &G, true);
        } else {
                LouvainMethod{}.performClustering(copy, &G, true);
        }

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = G.getPartitionIndex(node);
        } endfor

        G.set_partition_count(G.get_partition_count_compute());
        ind.objective     = ModularityMetric::computeModularity(G);
        ind.partition_map = partition_map;
        ind.cut_edges     = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor

        //if(output) {
                //m_filebuffer_string <<  m_global_timer.elapsed() <<  " " <<  ind.objective <<  std::endl;
                //m_time_stamp++;
        //}
}

void population_clustering::insert(graph_access & G, Individuum & ind) {
        if( ind.objective > best_objective ) {
                m_filebuffer_string <<  global_timer_elapsed() <<  " " <<  ind.objective <<  std::endl;
                m_time_stamp++;
                best_objective = ind.objective;
        }

        m_no_partition_calls++;
        if(m_internal_population_clustering.size() < m_population_clustering_size) {
                m_internal_population_clustering.push_back(ind);
        } else {
                double worst_objective = 1;
                for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                        if(m_internal_population_clustering[i].objective < worst_objective) {
                                worst_objective = m_internal_population_clustering[i].objective;
                        }
                }         
                if(ind.objective < worst_objective ) {
                        delete[] (ind.partition_map);
                        delete ind.cut_edges;
                        return; // do nothing
                }
                //else measure similarity
                unsigned max_similarity = std::numeric_limits<unsigned>::max();
                unsigned max_similarity_idx = 0;
                for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                        if(m_internal_population_clustering[i].objective <= ind.objective) {
                                //now measure
                                int diff_size = m_internal_population_clustering[i].cut_edges->size() + ind.cut_edges->size();
                                std::vector<EdgeID> output_diff(diff_size,std::numeric_limits<EdgeID>::max());

                                set_symmetric_difference(m_internal_population_clustering[i].cut_edges->begin(),
                                                m_internal_population_clustering[i].cut_edges->end(),
                                                ind.cut_edges->begin(),
                                                ind.cut_edges->end(),
                                                output_diff.begin());

                                unsigned similarity = 0;
                                for( unsigned j = 0; j < output_diff.size(); j++) {
                                        if(output_diff[j] < std::numeric_limits<EdgeID>::max()) {
                                                similarity++;
                                        } else {
                                                break;
                                        }
                                }

                                if( similarity < max_similarity) {
                                        max_similarity     = similarity;
                                        max_similarity_idx = i;
                                }
                        }
                }         

                delete[] (m_internal_population_clustering[max_similarity_idx].partition_map);
                delete m_internal_population_clustering[max_similarity_idx].cut_edges;

                m_internal_population_clustering[max_similarity_idx] = ind;
        }
}

void population_clustering::replace(Individuum & in, Individuum & out) {
        //first find it:
        for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                if(m_internal_population_clustering[i].partition_map == in.partition_map) {
                        //found it
                        delete[] (m_internal_population_clustering[i].partition_map);
                        delete m_internal_population_clustering[i].cut_edges;

                        m_internal_population_clustering[i] = out;
                        break;
                }
        }
}

void population_clustering::combine_basic_flat(const PartitionConfig & partition_config, 
                graph_access & G, 
                Individuum & first_ind, 
                Individuum & second_ind, 
                Individuum & output_ind) {

        std::vector< unsigned > lhs(G.number_of_nodes(), 0);
        std::vector< unsigned > rhs(G.number_of_nodes(), 0);
        std::vector< unsigned > output(G.number_of_nodes(), 0);

        forall_nodes(G, node) {
                lhs[node] = first_ind.partition_map[node];
                rhs[node] = second_ind.partition_map[node];
        } endfor

        output = maxmimum_overlap( lhs, rhs);
        graph_access contracted_graph = contract_by_clustering(G, output);

        clustering_t contracted_clustering; double quality;
        std::tie(contracted_clustering, quality) = do_louvain(contracted_graph);
        update_clustering(output, contracted_clustering);

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = output[node];
                G.setPartitionIndex(node, output[node]);
        } endfor

        G.set_partition_count(G.get_partition_count_compute());
        if (m_use_leiden) {
                splitDisconnectedCommunities(G);
                forall_nodes(G, n) { partition_map[n] = G.getPartitionIndex(n); } endfor
        }
        output_ind.objective = ModularityMetric::computeModularity(G);
        output_ind.partition_map    = partition_map;
        output_ind.cut_edges        = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                output_ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor
}

void population_clustering::combine_improved_multilevel(const PartitionConfig & partition_config,
                graph_access & G,
                Individuum & first_ind,
                Individuum & second_ind,
                Individuum & output_ind) {

        double eps = 0.001;
        std::vector< unsigned > lhs(G.number_of_nodes(), 0);
        std::vector< unsigned > rhs(G.number_of_nodes(), 0);
        std::vector< unsigned > output(G.number_of_nodes(), 0);

        forall_nodes(G, node) {
                lhs[node] = first_ind.partition_map[node];
                rhs[node] = second_ind.partition_map[node];
        } endfor

        output = maxmimum_overlap( lhs, rhs);

        graph_access* Q = new graph_access();
        G.copy(*Q);

        graph_hierarchy hierarchy;
        std::list<graph_access*> junk{ Q };

        std::vector< unsigned > current_clustering;
        std::vector< unsigned > contracted_overlap = output;

        //// constrained louvain
        double q = -1, q_;

        std::mt19937 gen{ std::random_device{}() };
        while(true) {
                current_clustering.resize(Q->number_of_nodes());
                std::iota(current_clustering.begin(), current_clustering.end(), 0);
                set_second_partition_index(*Q, contracted_overlap);

                q_ = q;
                q = local_search(*Q, current_clustering, gen, true);
                if(!(q - q_ > eps)) { break; }

                {
                        PartitionConfig config; // never read!
                        Q = Coarsening::performCoarsening(config, *Q, hierarchy, junk);
                }

                contracted_overlap = apply_fine_clustering_to_coarse_graph(contracted_overlap, current_clustering, Q->number_of_nodes());
        }

        clustering_t contracted_clustering;
        if( first_ind.objective > second_ind.objective ) {
                contracted_clustering = contract_better_clustering_by_contracted_overlap(output, contracted_overlap, lhs);
        } else {
                contracted_clustering = contract_better_clustering_by_contracted_overlap(output, contracted_overlap, rhs);
        }
        forall_nodes((*Q), node) {
                Q->setPartitionIndex(node, contracted_clustering[node]);
        } endfor
        Q->set_partition_count(Q->get_partition_count_compute());

        // apply the contracted clustering to the graph and push it onto the hierarchy
        hierarchy.push_back(Q, 0);

        // uncoarsen
        while(!hierarchy.isEmpty()) {
                Q = hierarchy.pop_finer_and_project();
                extract_clustering(*Q, current_clustering);

                q = local_search(*Q, current_clustering, gen);
        }

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = current_clustering[node];
                G.setPartitionIndex(node, current_clustering[node]);
        } endfor

        G.set_partition_count(G.get_partition_count_compute());
        if (m_use_leiden) {
                splitDisconnectedCommunities(G);
                forall_nodes(G, n) { partition_map[n] = G.getPartitionIndex(n); } endfor
        }
        output_ind.objective = ModularityMetric::computeModularity(G);
        output_ind.partition_map    = partition_map;
        output_ind.cut_edges        = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                output_ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor


        for(auto* g: junk) { delete g; }
}

void population_clustering::combine_improved_flat(const PartitionConfig & partition_config, 
                graph_access & G, 
                Individuum & first_ind, 
                Individuum & second_ind, 
                Individuum & output_ind) {

        std::vector< unsigned > lhs(G.number_of_nodes(), 0);
        std::vector< unsigned > rhs(G.number_of_nodes(), 0);
        std::vector< unsigned > output(G.number_of_nodes(), 0);

        forall_nodes(G, node) {
                lhs[node] = first_ind.partition_map[node];
                rhs[node] = second_ind.partition_map[node];
        } endfor

        output = maxmimum_overlap( lhs, rhs);
        graph_access contracted_graph = contract_by_clustering(G, output);

        clustering_t contracted_clustering;
        if( first_ind.objective > second_ind.objective ) {
                contracted_clustering = apply_fine_clustering_to_coarse_graph(lhs, output, contracted_graph.number_of_nodes());
        } else {
                contracted_clustering = apply_fine_clustering_to_coarse_graph(rhs, output, contracted_graph.number_of_nodes());
        }
        forall_nodes(contracted_graph, node) {
                contracted_graph.setPartitionIndex(node, contracted_clustering[node]);
        } endfor

        clustering_t new_contracted_clustering; double quality;
        std::tie(new_contracted_clustering, quality) = do_louvain(contracted_graph, contracted_clustering);
        update_clustering(output, new_contracted_clustering);

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = output[node];
                G.setPartitionIndex(node, output[node]);
        } endfor

        G.set_partition_count(G.get_partition_count_compute());
        if (m_use_leiden) {
                splitDisconnectedCommunities(G);
                forall_nodes(G, n) { partition_map[n] = G.getPartitionIndex(n); } endfor
        }
        output_ind.objective = ModularityMetric::computeModularity(G);
        output_ind.partition_map    = partition_map;
        output_ind.cut_edges        = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                output_ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor

}

void population_clustering::combine_improved_flat_with_partitioning(const PartitionConfig & partition_config, 
                graph_access & G, 
                Individuum & first_ind, 
                Individuum & output_ind) {


        std::vector< unsigned > lhs(G.number_of_nodes(), 0);
        std::vector< unsigned > rhs(G.number_of_nodes(), 0);
        std::vector< unsigned > output(G.number_of_nodes(), 0);

        forall_nodes(G, node) {
                lhs[node] = first_ind.partition_map[node];
        } endfor


        PartitionConfig cross_config;
        configuration{}.standard(cross_config);
        configuration{}.fastsocial(cross_config);

        int kfactor    = random_functions::nextInt(2,64);
        unsigned larger_imbalance = random_functions::nextInt(3,50);
        double epsilon = larger_imbalance/100.0;
        cross_config.k                                    = kfactor;
        cross_config.kaffpa_perfectly_balanced_refinement = false;
        cross_config.upper_bound_partition                = (1+epsilon)*ceil(G.number_of_nodes()/(double)kfactor);
        cross_config.refinement_scheduling_algorithm      = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
        cross_config.combine                              = false;
        cross_config.graph_allready_partitioned           = false;

        //G.set_partition_count(kfactor);
        balance_configuration bc;
        bc.configurate_balance( cross_config, G);

        graph_partitioner partitioner;
        partitioner.perform_partitioning(cross_config, G);

        forall_nodes(G, node) {
                rhs[node] = G.getPartitionIndex(node);
        } endfor

        output = maxmimum_overlap( lhs, rhs);
        graph_access contracted_graph = contract_by_clustering(G, output);

        clustering_t contracted_clustering = apply_fine_clustering_to_coarse_graph(lhs, output, contracted_graph.number_of_nodes());
        forall_nodes(contracted_graph, node) {
                contracted_graph.setPartitionIndex(node, contracted_clustering[node]);
        } endfor

        clustering_t new_contracted_clustering; double quality;
        std::tie(new_contracted_clustering, quality) = do_louvain(contracted_graph, contracted_clustering);
        update_clustering(output, new_contracted_clustering);

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = output[node];
                G.setPartitionIndex(node, output[node]);
        } endfor
        G.set_partition_count(G.get_partition_count_compute());

        output_ind.objective = ModularityMetric::computeModularity(G);
        output_ind.partition_map    = partition_map;
        output_ind.cut_edges        = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                output_ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor
}


void population_clustering::combine_improved_flat_with_sclp(const PartitionConfig & partition_config, 
                graph_access & G, 
                Individuum & first_ind, 
                Individuum & output_ind) {


        std::vector< unsigned > lhs(G.number_of_nodes(), 0);
        std::vector< unsigned > rhs(G.number_of_nodes(), 0);
        std::vector< unsigned > output(G.number_of_nodes(), 0);

        forall_nodes(G, node) {
                lhs[node] = first_ind.partition_map[node];
        } endfor


        PartitionConfig misc;
        configuration cfg;
        cfg.fastsocial(misc);

        misc.graph_allready_partitioned = false;
        misc.cluster_upperbound = random_functions::nextInt(10, G.number_of_nodes());
        misc.upper_bound_partition = misc.cluster_upperbound; 
        misc.cluster_coarsening_factor = 1;

        unsigned no_blocks = 0;
        size_constraint_label_propagation{}.label_propagation( misc, G, rhs, no_blocks);

        output = maxmimum_overlap( lhs, rhs);
        graph_access contracted_graph = contract_by_clustering(G, output);

        clustering_t contracted_clustering = apply_fine_clustering_to_coarse_graph(lhs, output, contracted_graph.number_of_nodes());
        forall_nodes(contracted_graph, node) {
                contracted_graph.setPartitionIndex(node, contracted_clustering[node]);
        } endfor

        clustering_t new_contracted_clustering; double quality;
        std::tie(new_contracted_clustering, quality) = do_louvain(contracted_graph, contracted_clustering);
        update_clustering(output, new_contracted_clustering);

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = output[node];
                G.setPartitionIndex(node, output[node]);
        } endfor

        G.set_partition_count(G.get_partition_count_compute());
        if (m_use_leiden) {
                splitDisconnectedCommunities(G);
                forall_nodes(G, n) { partition_map[n] = G.getPartitionIndex(n); } endfor
        }
        output_ind.objective = ModularityMetric::computeModularity(G);
        output_ind.partition_map    = partition_map;
        output_ind.cut_edges        = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                output_ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor
}

void population_clustering::mutate( const PartitionConfig & partition_config, graph_access & G, Individuum & first_ind, Individuum & second_ind, Individuum & output_ind) {
        Individuum output_a;
        Individuum output_b;
        mutate_random(partition_config, G, first_ind, output_a);
        mutate_random(partition_config, G, second_ind, output_b);
        combine_improved_multilevel(partition_config, G, output_a, output_b, output_ind);

        delete[] (output_a.partition_map);
        delete output_a.cut_edges;

        delete[] (output_b.partition_map);
        delete output_b.cut_edges;
}

void population_clustering::mutate_random( const PartitionConfig & partition_config, graph_access & G, Individuum & first_ind, Individuum & output_ind) {
        std::vector< unsigned > clustering(G.number_of_nodes(), 0);
        std::mt19937 gen{ std::random_device{}() };
        forall_nodes(G, node) {
                clustering[node] = first_ind.partition_map[node];
        } endfor

        double l = random_functions::nextDouble(0.01,partition_config.mh_mutate_fraction);
        unsigned c = *std::max_element(clustering.begin(), clustering.end()) + 1;
        unsigned clusters_to_select = std::ceil(l * c);

        //std::cout << "mutating with l = " << l << " (clusters_to_select = " << clusters_to_select << ")\n";

        assert(c != 0);
        std::uniform_int_distribution<size_t> dist{ 0, c - 1 };
        std::set<unsigned> selected_clusters;
        for(size_t i = 0; i < clusters_to_select; ++i) {
                while(!selected_clusters.insert(dist(gen)).second);
        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, clustering[node]);
        } endfor
        G.set_partition_count(G.get_partition_count_compute());

        std::uniform_real_distribution<double> dist_eps{ 0.1, 0.5 };
        for(unsigned cluster: selected_clusters) {
                graph_access E;
                std::vector<unsigned> mapping;
                graph_extractor{ }.extract_block(G, E, cluster, mapping);

                PartitionConfig partition_config;

                configuration cfg;
                partition_config.k         = 2;
                cfg.standard(partition_config);

                double eps                 = dist_eps(gen);
                partition_config.imbalance = 100*eps;
                partition_config.epsilon   = 100*eps;

                cfg.fast(partition_config);

                E.set_partition_count(2);
                balance_configuration bc;
                bc.configurate_balance(partition_config, E);

                graph_partitioner{ }.perform_partitioning(partition_config, E);
                for(size_t i = 0; i < mapping.size(); ++i) {
                        if(E.getPartitionIndex(i) == 1) { clustering[mapping[i]] = c; }
                }
                ++c;
        }

        local_search(G, clustering, gen);

        int* partition_map = new int[G.number_of_nodes()];
        forall_nodes(G, node) {
                partition_map[node] = clustering[node];
                G.setPartitionIndex(node, clustering[node]);
        } endfor

        G.set_partition_count(G.get_partition_count_compute());
        if (m_use_leiden) {
                splitDisconnectedCommunities(G);
                forall_nodes(G, n) { partition_map[n] = G.getPartitionIndex(n); } endfor
        }
        output_ind.objective = ModularityMetric::computeModularity(G);
        output_ind.partition_map    = partition_map;
        output_ind.cut_edges        = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                output_ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor


}

void population_clustering::extinction( ) {
        for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                delete[] m_internal_population_clustering[i].partition_map;
                delete m_internal_population_clustering[i].cut_edges; 
        }

        m_internal_population_clustering.clear();
        m_internal_population_clustering.resize(0);
}

void population_clustering::get_two_random_individuals(Individuum & first, Individuum & second) {
        int first_idx = random_functions::nextInt(0, m_internal_population_clustering.size()-1);
        first = m_internal_population_clustering[first_idx];

        int second_idx = random_functions::nextInt(0, m_internal_population_clustering.size()-1);
        while( first_idx == second_idx ) {
                second_idx = random_functions::nextInt(0, m_internal_population_clustering.size()-1);
        }

        second = m_internal_population_clustering[second_idx];
}

void population_clustering::get_one_individual_tournament(Individuum & first) {
        Individuum one, two;
        get_two_random_individuals(one, two);
        first  =  one.objective > two.objective ? one : two;
}

void population_clustering::get_two_individuals_tournament(Individuum & first, Individuum & second) {
        Individuum one, two;
        get_two_random_individuals(one, two);
        first  =  one.objective > two.objective? one : two;

        get_two_random_individuals(one, two);
        second =  one.objective > two.objective ? one : two;

        if( first.objective == second.objective) {
                second = one.objective >= two.objective? one : two;
        }
}

void population_clustering::get_random_individuum(Individuum & ind) {
        int idx = random_functions::nextInt(0, m_internal_population_clustering.size()-1);
        ind     = m_internal_population_clustering[idx];
}

void population_clustering::get_best_individuum(Individuum & ind) {
        double max_objective = -1;
        unsigned idx         = 0;

        for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                if(m_internal_population_clustering[i].objective > max_objective) {
                        max_objective = m_internal_population_clustering[i].objective;
                        idx           = i;
                }
        }

        ind = m_internal_population_clustering[idx];
}

bool population_clustering::is_full() {
        return m_internal_population_clustering.size() == m_population_clustering_size;
}

void population_clustering::apply_fittest( graph_access & G, double & objective ) {
        double max_objective = -1;
        unsigned idx         = 0;

        quality_metrics qm;
        for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, m_internal_population_clustering[i].partition_map[node]);
                } endfor
                if(m_internal_population_clustering[i].objective > max_objective) {
                        max_objective = m_internal_population_clustering[i].objective;
                        idx           = i;
                }
        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, m_internal_population_clustering[idx].partition_map[node]);
        } endfor

        objective = max_objective;
}

void population_clustering::print() {
        int rank;
        MPI_Comm_rank( m_communicator, &rank);

        std::cout <<  "rank " <<  rank << " fingerprint ";

        for( unsigned i = 0; i < m_internal_population_clustering.size(); i++) {
                std::cout <<  m_internal_population_clustering[i].objective << " ";
        }         

        std::cout <<  std::endl;
}

void population_clustering::write_log(std::string & filename) {
        m_filebuffer_string <<  global_timer_elapsed() <<  " " <<  best_objective <<  std::endl;
        std::ofstream f(filename.c_str());
        f << m_filebuffer_string.str();
        f.close();
}

