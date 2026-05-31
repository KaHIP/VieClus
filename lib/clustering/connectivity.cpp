/******************************************************************************
 * connectivity.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/


#include "connectivity.h"

#include <queue>
#include <unordered_map>
#include <vector>

PartitionID splitDisconnectedCommunities(graph_access &G) {
        NodeID n = G.number_of_nodes();
        PartitionID num_clusters = G.get_partition_count();

        // build per-cluster node lists
        std::vector<std::vector<NodeID>> cluster_nodes(num_clusters);
        forall_nodes(G, node) {
                cluster_nodes[G.getPartitionIndex(node)].push_back(node);
        } endfor

        std::vector<bool> visited(n, false);
        PartitionID next_id = num_clusters;

        for (PartitionID c = 0; c < num_clusters; ++c) {
                if (cluster_nodes[c].empty()) continue;

                unsigned component_count = 0;
                for (NodeID seed : cluster_nodes[c]) {
                        if (visited[seed]) continue;

                        // BFS within cluster c
                        std::queue<NodeID> bfs_queue;
                        std::vector<NodeID> component;
                        bfs_queue.push(seed);
                        visited[seed] = true;

                        while (!bfs_queue.empty()) {
                                NodeID cur = bfs_queue.front();
                                bfs_queue.pop();
                                component.push_back(cur);

                                forall_out_edges(G, e, cur) {
                                        NodeID target = G.getEdgeTarget(e);
                                        if (!visited[target] && G.getPartitionIndex(target) == c) {
                                                visited[target] = true;
                                                bfs_queue.push(target);
                                        }
                                } endfor
                        }

                        // first component keeps original ID, others get new IDs
                        if (component_count > 0) {
                                for (NodeID v : component) {
                                        G.setPartitionIndex(v, next_id);
                                }
                                next_id++;
                        }
                        component_count++;
                }
        }

        // remap to consecutive IDs [0, k-1]
        std::unordered_map<PartitionID, PartitionID> new_mapping;
        PartitionID id = 0;

        forall_nodes(G, node) {
                PartitionID cluster = G.getPartitionIndex(node);
                if (!new_mapping.count(cluster)) { new_mapping[cluster] = id++; }
                G.setPartitionIndex(node, new_mapping[cluster]);
        } endfor

        G.set_partition_count(id);
        return id;
}
