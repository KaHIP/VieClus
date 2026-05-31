/******************************************************************************
 * leidenmethod.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/


#include "leidenmethod.h"

#include "clustering/coarsening/coarsening.h"
#include "clustering/labelpropagation.h"
#include "clustering/neighborhood.h"
#include "clustering/connectivity.h"
#include "coarsening/clustering/size_constraint_label_propagation.h"
#include "partition/coarsening/clustering/node_ordering.h"
#include "tools/modularitymetric.h"
#include "tools/random_functions.h"

#include <algorithm>
#include <cmath>
#include <list>
#include <numeric>
#include <unordered_map>
#include <vector>

using namespace std;

LeidenMethod::LeidenMethod()
    : m_G(0)
{
}

LeidenMethod::~LeidenMethod()
{
}

PartitionID LeidenMethod::performClustering(PartitionConfig &config, graph_access *G, bool start_w_singletons)
{
    return LeidenMethod::performClusteringWithLPP(config, G, start_w_singletons);
}

PartitionID LeidenMethod::performClusteringWithLPP(const PartitionConfig& config,
                                                    graph_access* G, bool start_w_singletons)
{
    NodeID numberOfMoves = 0;
    unsigned coarsenings = 0;
    graph_hierarchy graphHierarchy;
    list<graph_access *> coarseGraphsToDelete;

    m_G = G;

    // optional label propagation preprocessing (same as Louvain)
    for (unsigned i = 0; i < config.lm_number_of_label_propagation_levels; ++i)
    {
        if (start_w_singletons) { initializeSingletonClusters(); }

        if (config.lm_cluster_coarsening_factor > 0)
        {
            NodeID no_blocks = 0;
            std::vector<NodeID> cluster_id;
            size_constraint_label_propagation sclp;
            sclp.label_propagation(config, *m_G, cluster_id, no_blocks);

            forall_nodes((*m_G), node) {
                m_G->setPartitionIndex(node, cluster_id[node]);
            } endfor

            // enforce connectivity after size-constrained LP
            splitDisconnectedCommunities(*m_G);

            numberOfMoves = 1;
        }
        else
        {
            numberOfMoves = LabelPropagation::performLabelPropagation(config, m_G);
            // enforce connectivity after LP
            splitDisconnectedCommunities(*m_G);
        }

        if (numberOfMoves)
        {
            m_G = Coarsening::performCoarsening(config, *m_G, graphHierarchy, coarseGraphsToDelete);
            coarsenings++;
        }

        if (!numberOfMoves)
        {
            break;
        }
    }

    // main Leiden loop: Phase 1 (fast local move) + Phase 2 (refinement) + Phase 3 (aggregation)
    do
    {
        if (start_w_singletons) { initializeSingletonClusters(); }

        // Phase 1: fast local moving (queue-based)
        numberOfMoves = performFastLocalMove(config);

        if (numberOfMoves)
        {
            // save the flat (non-refined) partition from Phase 1
            vector<PartitionID> flatPartition(m_G->number_of_nodes());
            forall_nodes((*m_G), node) {
                flatPartition[node] = m_G->getPartitionIndex(node);
            } endfor

            // Phase 2: refinement (guarantees connected communities)
            performRefinement(config, flatPartition);

            // Phase 3: aggregate based on refined partition
            m_G = Coarsening::performCoarsening(config, *m_G, graphHierarchy, coarseGraphsToDelete);
            coarsenings++;

            // Initialize next level with flat (non-refined) partition projected onto coarse graph.
            // The coarse graph nodes correspond to refined sub-communities.
            // We need to map each coarse node back to its flat partition community.
            // After coarsening, partition indices on coarse graph are initialized by the contractor.
            // We override them with the flat partition mapping.
            //
            // Each node in the coarse graph represents a refined sub-community.
            // We find which flat community it belongs to by looking at the mapping.
            CoarseMapping* mapping = graphHierarchy.get_mapping_of_current_finer();
            if (mapping != nullptr) {
                // mapping[fine_node] = coarse_node
                // flatPartition[fine_node] = flat community
                // For each coarse node, pick the flat community of any constituent fine node
                NodeID coarse_n = m_G->number_of_nodes();
                vector<PartitionID> coarse_flat(coarse_n, 0);
                vector<bool> assigned(coarse_n, false);

                NodeID fine_n = flatPartition.size();
                for (NodeID fine_node = 0; fine_node < fine_n; ++fine_node) {
                    NodeID coarse_node = (*mapping)[fine_node];
                    if (!assigned[coarse_node]) {
                        coarse_flat[coarse_node] = flatPartition[fine_node];
                        assigned[coarse_node] = true;
                    }
                }

                // remap flat community IDs to consecutive range
                unordered_map<PartitionID, PartitionID> flat_remap;
                PartitionID fid = 0;
                for (NodeID cn = 0; cn < coarse_n; ++cn) {
                    if (!flat_remap.count(coarse_flat[cn])) {
                        flat_remap[coarse_flat[cn]] = fid++;
                    }
                    m_G->setPartitionIndex(cn, flat_remap[coarse_flat[cn]]);
                }
                m_G->set_partition_count(fid);

                // next iteration will NOT start with singletons since we initialized
                // from the flat partition
                start_w_singletons = false;
            }
        }
    }
    while (numberOfMoves);

    // append the last level
    graphHierarchy.push_back(m_G, 0);

    // uncoarsening with refinement
    while (coarsenings > 0 && !graphHierarchy.isEmpty())
    {
        m_G = graphHierarchy.pop_finer_and_project();
        numberOfMoves = performFastLocalMove(config);
    }

    // cleanup coarse graphs
    while (!coarseGraphsToDelete.empty())
    {
        delete coarseGraphsToDelete.front();
        coarseGraphsToDelete.pop_front();
    }

    // final connectivity enforcement
    splitDisconnectedCommunities(*m_G);

    // remap cluster IDs to consecutive range
    std::unordered_map<PartitionID, PartitionID> new_mapping;
    PartitionID id = 0;

    forall_nodes((*m_G), node) {
        PartitionID cluster = m_G->getPartitionIndex(node);
        if (!new_mapping.count(cluster)) { new_mapping[cluster] = id++; }
        m_G->setPartitionIndex(node, new_mapping[cluster]);
    } endfor

    m_G->set_partition_count(id);

    return m_G->get_partition_count();
}


void LeidenMethod::initializeSingletonClusters()
{
    PartitionID clusterID = 0;
    forall_nodes((*m_G), n)
    {
        m_G->setPartitionIndex(n, clusterID);
        clusterID++;
    } endfor
    m_G->set_partition_count(clusterID);
}


NodeID LeidenMethod::performFastLocalMove(const PartitionConfig &config)
{
    NodeID numberOfMoves = 0;
    bool hasGraphSelfLoops = m_G->containsSelfLoops();
    Neighborhood neighborhood;
    neighborhood.initialize(m_G);

    ModularityMetric *objective = new ModularityMetric(*m_G);

    NodeID n = m_G->number_of_nodes();

    // queue-based: start with all nodes in random order
    vector<NodeID> perm(n);
    random_functions::permutate_vector_good(perm, true);

    queue<NodeID> node_queue;
    vector<bool> in_queue(n, false);

    for (NodeID i = 0; i < n; ++i) {
        node_queue.push(perm[i]);
        in_queue[perm[i]] = true;
    }

    while (!node_queue.empty())
    {
        NodeID node = node_queue.front();
        node_queue.pop();
        in_queue[node] = false;

        neighborhood.update(node);

        if (neighborhood.getNumberOfNeighboringClusters() > 1)
        {
            PartitionID oldCluster = m_G->getPartitionIndex(node);
            PartitionID bestCluster = oldCluster;
            double bestGain = 0.0;
            EdgeWeight selfLoop = 0;

            if (hasGraphSelfLoops)
            {
                selfLoop = m_G->getSelfLoop(node);
            }

            objective->removeNode(node, oldCluster, neighborhood.getEdgeWeightToNeighboringCluster(oldCluster), selfLoop);

            for (NodeID i = 0; i < neighborhood.getNumberOfNeighboringClusters(); ++i)
            {
                PartitionID newCluster = neighborhood.getClusterIDOfNeighbor(i);
                double gain = objective->gain(node, newCluster, neighborhood.getEdgeWeightToNeighboringCluster(newCluster));

                if (bestGain < gain)
                {
                    bestGain = gain;
                    bestCluster = newCluster;
                }
            }

            objective->insertNode(node, bestCluster, neighborhood.getEdgeWeightToNeighboringCluster(bestCluster), selfLoop);

            if (oldCluster != bestCluster)
            {
                numberOfMoves++;
                // enqueue neighbors that are not already queued
                forall_out_edges((*m_G), e, node) {
                    NodeID neighbor = m_G->getEdgeTarget(e);
                    if (!in_queue[neighbor]) {
                        node_queue.push(neighbor);
                        in_queue[neighbor] = true;
                    }
                } endfor
            }
        }
    }

    delete objective;
    return numberOfMoves;
}


void LeidenMethod::performRefinement(const PartitionConfig &config,
                                      const std::vector<PartitionID> &flatPartition)
{
    NodeID n = m_G->number_of_nodes();
    bool hasGraphSelfLoops = m_G->containsSelfLoops();

    // build per-community node lists from the flat partition
    PartitionID numFlatClusters = *max_element(flatPartition.begin(), flatPartition.end()) + 1;
    vector<vector<NodeID>> communityNodes(numFlatClusters);
    for (NodeID node = 0; node < n; ++node) {
        communityNodes[flatPartition[node]].push_back(node);
    }

    // reset all nodes to singletons (refined partition starts from scratch)
    initializeSingletonClusters();

    // track which nodes are still singletons (have not been merged)
    vector<bool> isSingleton(n, true);

    // create modularity objective for the singleton state
    ModularityMetric *objective = new ModularityMetric(*m_G);

    double theta = m_leiden_config.theta;

    // for each community in the flat partition, refine
    for (PartitionID c = 0; c < numFlatClusters; ++c)
    {
        vector<NodeID>& nodes = communityNodes[c];
        if (nodes.size() <= 1) continue;

        // shuffle nodes within this community
        for (size_t i = nodes.size() - 1; i > 0; --i) {
            size_t j = random_functions::nextInt(0, i);
            swap(nodes[i], nodes[j]);
        }

        for (NodeID node : nodes)
        {
            if (!isSingleton[node]) continue;

            // find adjacent refined sub-communities within the same flat community
            PartitionID myCluster = m_G->getPartitionIndex(node);
            EdgeWeight selfLoop = hasGraphSelfLoops ? m_G->getSelfLoop(node) : 0;

            // collect candidates: neighboring refined clusters within same flat community C
            struct Candidate {
                PartitionID cluster;
                EdgeWeight edgeWeight;
                double deltaQ;
            };
            vector<Candidate> candidates;
            double totalWeight = 0.0;

            // compute edge weights to neighboring refined clusters within C
            unordered_map<PartitionID, EdgeWeight> edgeWeightToCluster;
            forall_out_edges((*m_G), e, node) {
                NodeID neighbor = m_G->getEdgeTarget(e);
                if (flatPartition[neighbor] == c) {
                    PartitionID neighborCluster = m_G->getPartitionIndex(neighbor);
                    if (neighborCluster != myCluster) {
                        edgeWeightToCluster[neighborCluster] += m_G->getEdgeWeight(e);
                    }
                }
            } endfor

            if (edgeWeightToCluster.empty()) continue;

            // temporarily remove node from its singleton cluster
            EdgeWeight edgeWeightToOwn = 0;
            forall_out_edges((*m_G), e, node) {
                NodeID neighbor = m_G->getEdgeTarget(e);
                if (m_G->getPartitionIndex(neighbor) == myCluster) {
                    edgeWeightToOwn += m_G->getEdgeWeight(e);
                }
            } endfor

            objective->removeNode(node, myCluster, edgeWeightToOwn, selfLoop);

            for (auto& kv : edgeWeightToCluster) {
                double gain = objective->gain(node, kv.first, kv.second);
                if (gain > 0) {
                    double weight = exp(gain / theta);
                    candidates.push_back({kv.first, kv.second, gain});
                    totalWeight += weight;
                }
            }

            if (candidates.empty()) {
                // put node back in its own cluster
                objective->insertNode(node, myCluster, edgeWeightToOwn, selfLoop);
                continue;
            }

            // stochastic selection proportional to exp(deltaQ / theta)
            double r = random_functions::nextDouble(0.0, totalWeight);
            double cumulative = 0.0;
            PartitionID selectedCluster = candidates[0].cluster;
            EdgeWeight selectedEdgeWeight = candidates[0].edgeWeight;

            for (auto& cand : candidates) {
                cumulative += exp(cand.deltaQ / theta);
                if (cumulative >= r) {
                    selectedCluster = cand.cluster;
                    selectedEdgeWeight = cand.edgeWeight;
                    break;
                }
            }

            // merge node into selected sub-community
            objective->insertNode(node, selectedCluster, selectedEdgeWeight, selfLoop);
            isSingleton[node] = false;
        }
    }

    delete objective;
}
