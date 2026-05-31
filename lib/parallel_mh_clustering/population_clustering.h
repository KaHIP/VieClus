/******************************************************************************
 * population_clustering.h
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/

#ifndef POPULATION_AEFH46G6A
#define POPULATION_AEFH46G6A

#include <sstream>
#ifdef USE_MPI
#include <mpi.h>
#else
#include "tools/pseudo_mpi.h"
#endif
#include <unordered_map>
#include <functional>
#include <algorithm>

#include "data_structure/graph_access.h"
#include "clustering/coarsening/coarsening.h"
#include "partition_config.h"
#include "tools/global_timer.h"
#include "clustering/louvainmethod.h"
#include "clustering/leiden_config.h"
#include "configuration.h"
#include "tools/modularitymetric.h"
#include "tools/random_functions.h"

struct Individuum {
        int* partition_map;
        double objective;
        std::vector<EdgeID>* cut_edges; //sorted
};

struct ENC {
        std::vector<NodeID> vertices;
};

using clustering_t = std::vector<unsigned>;

class population_clustering {
        public:
                population_clustering( MPI_Comm comm, const PartitionConfig & config, const LeidenConfig & leiden_config = LeidenConfig() );
                virtual ~population_clustering();

                void createIndividuum(const PartitionConfig & config, 
                                graph_access & G, 
                                Individuum & ind, 
                                bool output); 

                void combine_basic_flat(const PartitionConfig & config, 
                                graph_access & G, 
                                Individuum & first_ind, 
                                Individuum & second_ind, 
                                Individuum & output_ind); 

                void combine_improved_multilevel(const PartitionConfig & config, 
                                graph_access & G, 
                                Individuum & first_ind, 
                                Individuum & second_ind, 
                                Individuum & output_ind); 

                void combine_improved_flat(const PartitionConfig & config, 
                                graph_access & G, 
                                Individuum & first_ind, 
                                Individuum & second_ind, 
                                Individuum & output_ind); 

                void combine_improved_flat_with_sclp(const PartitionConfig & config, 
                                graph_access & G, 
                                Individuum & first_ind, 
                                Individuum & output_ind); 

                void combine_improved_flat_with_partitioning(const PartitionConfig & config, 
                                graph_access & G, 
                                Individuum & first_ind, 
                                Individuum & output_ind); 

                void combine_cross(const PartitionConfig & partition_config, 
                                graph_access & G, 
                                Individuum & first_ind, 
                                Individuum & output_ind);

                void mutate_random(const PartitionConfig & partition_config, 
                                graph_access & G, 
                                Individuum & first_ind,
                                Individuum & output_ind);

                void mutate( const PartitionConfig & partition_config, 
                             graph_access & G, 
                             Individuum & first_ind, 
                             Individuum & second_ind, 
                             Individuum & output_ind);

                void insert(graph_access & G, Individuum & ind);

                void set_pool_size(int size);

                void extinction();

                void get_two_random_individuals(Individuum & first, Individuum & second);

                void get_one_individual_tournament(Individuum & first); 

                void get_two_individuals_tournament(Individuum & first, Individuum & second);

                void replace(Individuum & in, Individuum & out);

                void get_random_individuum(Individuum & ind);

                void get_best_individuum(Individuum & ind);

                bool is_full(); 

                void apply_fittest( graph_access & G, double & objective);

                unsigned size() { return m_internal_population_clustering.size(); }

                void print();

                void write_log(std::string & filename);

                /* updates old clustering with a new clustering of a possibly contracted graph.
                 * new_coarse_clustering must be of size max(clustering) + 1, i.e. must have an
                 * entry for every cluster in the old clustering. */
                void update_clustering(clustering_t& clustering, clustering_t& new_coarse_clustering) {
                        canonicalize(new_coarse_clustering);
                        for(auto& cluster: clustering) { cluster = new_coarse_clustering[cluster]; }
                }

                clustering_t maxmimum_overlap(clustering_t const& a, clustering_t const& b) {
                        assert(a.size() == b.size() && "a and b of unequal length");

                        using ensemble_tuple = std::tuple<size_t, unsigned, unsigned>;
                        std::unordered_map<ensemble_tuple, unsigned, std::function<size_t(ensemble_tuple)>, std::function<bool(ensemble_tuple, ensemble_tuple)>> mapping(
                                        10,
                                        [](ensemble_tuple const& t) { return std::get<0>(t) * std::get<1>(t) + std::get<2>(t); },
                                        [](ensemble_tuple const& t, ensemble_tuple const& u) { return std::get<1>(t) == std::get<1>(u) && std::get<2>(t) == std::get<2>(u); }
                                        );

                        clustering_t overlap(a.size());
                        size_t cluster_id = 0;
                        size_t const n = a.size();

                        for(unsigned vertex = 0; vertex < n; ++vertex) {
                                auto tup = std::make_tuple(n, a[vertex], b[vertex]);
                                if(!mapping.count(tup)) { mapping[tup] = cluster_id++; }
                                overlap[vertex] = mapping[tup];
                        }

                        return overlap;
                }

                graph_access contract_by_clustering(graph_access& G, clustering_t const& clustering) {
                        graph_access G_;
                        G.copy(G_);

                        apply_clustering(G_, clustering);

                        PartitionConfig config; // never read!
                        graph_hierarchy gh;
                        std::list<graph_access*> junk;

                        graph_access* G_coarse = Coarsening::performCoarsening(config, G_, gh, junk);
                        G_coarse->copy(G_);
                        delete G_coarse;

                        return G_;
                }

                void apply_clustering(graph_access& G, clustering_t const& clustering) {
                        assert(G.number_of_nodes() == clustering.size());
                        forall_nodes(G, vertex) {
                                G.setPartitionIndex(vertex, clustering[vertex]);
                        } endfor

                        PartitionID clusters_num = *std::max_element(clustering.begin(), clustering.end()) + 1;
                        G.set_partition_count(clusters_num);
                }

                /* extracts the clustering stored in the partition indices of all nodes into a
                 * vector. */
                void extract_clustering(graph_access& G, clustering_t& clustering) {
                        clustering.resize(G.number_of_nodes());
                        forall_nodes(G, vertex) {
                                clustering[vertex] = G.getPartitionIndex(vertex);
                        } endfor

                        canonicalize(clustering);
                }

                void extract_second_clustering(graph_access& G, clustering_t& clustering) {
                        clustering.resize(G.number_of_nodes());
                        forall_nodes(G, vertex) {
                                clustering[vertex] = G.getSecondPartitionIndex(vertex);
                        } endfor

                        canonicalize(clustering);
                }

                /* remaps cluster ids to a canonical order */
                void canonicalize(clustering_t& clustering)  {
                        std::unordered_map<int, int> new_mapping;
                        int id = 0;

                        for(auto& cluster: clustering) {
                                if(new_mapping.count(cluster)) {
                                        cluster = new_mapping[cluster];
                                } else {
                                        cluster = new_mapping[cluster] = id++;
                                }
                        }
                }

                /* executes the louvain algorithm on the given graph in all its multilevel
                 * glory, and returns the found clustering along with its found quality. */
                std::pair<clustering_t, double> do_louvain(graph_access& G, clustering_t const& c = clustering_t{}) {
                        assert((!c.size() || G.number_of_nodes() == c.size()) && "clustering <-> graph size mismatch");

                        //srand(rand());
                        PartitionConfig partition_config;
                        configuration{}.standard(partition_config);

                        partition_config.lm_number_of_label_propagation_levels = random_functions::nextInt(0,5); 

                        partition_config.upper_bound_partition = G.number_of_nodes() + 1;
                        partition_config.graph_allready_partitioned = false;
                        partition_config.node_ordering = RANDOM_NODEORDERING;
                        partition_config.cluster_upperbound = random_functions::nextInt(G.number_of_nodes()/10, G.number_of_nodes());
                        partition_config.upper_bound_partition = partition_config.cluster_upperbound; 
                        partition_config.cluster_coarsening_factor = 1;


                        if (m_use_leiden) {
                                do_leiden_clustering(partition_config, G, c.empty());
                        } else {
                                LouvainMethod{}.performClustering(partition_config, &G, c.empty());
                        }

                        clustering_t clustering(G.number_of_nodes(), -1);
                        extract_clustering(G, clustering);

                        return {clustering, ModularityMetric::computeModularity(G)};
                }


                /* applies a fine (presumably the better one of the two going into overlap, but
                 * doesn't matter) onto a contracted clustering resulting from overlap.  e.g.
                 * if fine_good = [0, 0, 1, 1, 0] and overlap = [0, 0, 1, 1, 2], then the
                 * resulting contracted clustering is [0, 1, 0] since cluster 2 is part of
                 * cluster 0 in fine_good. */
                clustering_t apply_fine_clustering_to_coarse_graph(clustering_t const& fine_good, clustering_t const& overlap, size_t c) {
                        assert(fine_good.size() == overlap.size() && "size mismatch");
                        clustering_t coarse_clustering(c);
                        std::iota(coarse_clustering.begin(), coarse_clustering.end(), 0);

                        for(size_t i = 0; i < fine_good.size(); ++i) {
                                if(fine_good[i] != overlap[i]) { coarse_clustering[overlap[i]] = fine_good[i]; }
                        }

                        return coarse_clustering;
                } 

                void set_second_partition_index(graph_access& G, std::vector< NodeID > const& clustering) {
                        assert(G.number_of_nodes() == clustering.size());
                        G.resizeSecondPartitionIndex(clustering.size());

                        forall_nodes(G, node) {
                                G.setSecondPartitionIndex(node, clustering[node]);
                        } endfor
                }

                template <class RNG>
                        double local_search(graph_access& G, clustering_t& clustering, RNG&& gen, bool combine = false, double eps = 0.0001) {
                                assert(G.number_of_nodes() == clustering.size() && "no!");

                                apply_clustering(G, clustering);
                                ModularityMetric mod{ G };

                                std::vector<size_t> order(G.number_of_nodes());
                                std::iota(order.begin(), order.end(), 0);
                                std::shuffle(order.begin(), order.end(), gen);

                                double q = mod.quality(), q_;
                                do {
                                        for(size_t i = 0; i < G.number_of_nodes(); ++i) {
                                                size_t vertex        = order[i];
                                                unsigned cur_cluster = G.getPartitionIndex(vertex);

                                                std::unordered_map<size_t, size_t> hood_edges{ {cur_cluster, 0} };
                                                forall_out_edges(G, e, vertex) {
                                                        NodeID neighbor = G.getEdgeTarget(e);
                                                        size_t neighbor_cluster = G.getPartitionIndex(neighbor);

                                                        if(combine && G.getSecondPartitionIndex(vertex) != G.getSecondPartitionIndex(neighbor)) { continue; }

                                                        if(hood_edges.count(neighbor_cluster)) {
                                                                hood_edges[neighbor_cluster] += G.getEdgeWeight(e);
                                                        } else {
                                                                hood_edges[neighbor_cluster] = G.getEdgeWeight(e);
                                                        }
                                                } endfor

                                                double best_increase = 0;
                                                auto best_candidate = hood_edges.find(cur_cluster);

                                                mod.removeNode(vertex, cur_cluster, hood_edges[cur_cluster]);
                                                for(auto it = hood_edges.begin(); it != hood_edges.end(); ++it) {
                                                        double new_increase = mod.gain(vertex, it->first, it->second);
                                                        if(new_increase > best_increase) { best_candidate = it; best_increase = new_increase; }
                                                }

                                                mod.insertNode(vertex, best_candidate->first, best_candidate->second);
                                        }

                                        q_ = q; q = mod.quality();
                                } while(q - q_ > eps);

                                extract_clustering(G, clustering);
                                canonicalize(clustering);

                                return q;
                        }

                clustering_t contract_better_clustering_by_contracted_overlap(clustering_t const& overlap, clustering_t const& contracted_overlap, clustering_t const& better) {
                        std::unordered_map<NodeID, NodeID> mapping;
                        for(size_t i = 0; i < better.size(); ++i) {
                                if(!mapping.count(overlap[i])) { mapping[overlap[i]] = better[i]; }
                                assert(mapping[overlap[i]] == better[i]);
                        }

                        clustering_t contracted_better(contracted_overlap.size());
                        for(size_t i = 0; i < contracted_overlap.size(); ++i) { contracted_better[i] = mapping[contracted_overlap[i]]; }
                        return contracted_better;
                }

                void do_leiden_clustering(PartitionConfig & config, graph_access & G, bool start_w_singletons);

        private:
                bool m_use_leiden;
                double m_leiden_theta;

                unsigned                m_no_partition_calls;
                unsigned 		m_population_clustering_size;
                std::vector<Individuum> m_internal_population_clustering;
                std::vector< std::vector< unsigned int > > m_vertex_ENCs;
                std::vector< ENC > m_ENCs;

                int m_num_NCs;
                int m_num_NCs_computed;
                int m_num_ENCs;
                int m_time_stamp;
                double best_objective;

                MPI_Comm m_communicator;

                std::stringstream m_filebuffer_string;
};


#endif /* end of include guard: POPULATION_AEFH46G6 */
