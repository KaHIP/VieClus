/******************************************************************************
 * vieclus.cpp -- pybind11 wrapper for VieClus
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/

#include <pybind11/pybind11.h>
#include "interface/vieclus_interface.h"

pybind11::object wrap_vieclus(
                const pybind11::object &vwgt,
                const pybind11::object &xadj,
                const pybind11::object &adjwgt,
                const pybind11::object &adjncy,
                bool suppress_output,
                int seed,
                double time_limit,
                int cluster_upperbound,
                bool leiden,
                double leiden_theta) {
        int n = pybind11::len(xadj) - 1;
        std::vector<int> xadjv, adjncyv, vwgtv, adjwgtv;

        for (auto it : xadj)
                xadjv.push_back(pybind11::cast<int>(*it));

        for (auto it : adjncy)
                adjncyv.push_back(pybind11::cast<int>(*it));

        for (auto it : vwgt)
                vwgtv.push_back(pybind11::cast<int>(*it));

        for (auto it : adjwgt)
                adjwgtv.push_back(pybind11::cast<int>(*it));

        int* clustering   = new int[n];
        double modularity = 0;
        int num_clusters  = 0;

        vieclus_clustering(&n, &vwgtv[0], &xadjv[0],
                           &adjwgtv[0], &adjncyv[0],
                           suppress_output, seed,
                           time_limit, cluster_upperbound,
                           leiden, leiden_theta,
                           &modularity, &num_clusters, clustering);

        pybind11::list clustering_list;
        for (int i = 0; i < n; ++i)
                clustering_list.append(clustering[i]);

        delete[] clustering;

        return pybind11::make_tuple(modularity, clustering_list);
}

PYBIND11_MODULE(vieclus, m) {
        m.doc() = "VieClus -- Vienna Graph Clustering";
        m.def("cluster", &wrap_vieclus, "Cluster a graph using VieClus.",
              pybind11::arg("vwgt"),
              pybind11::arg("xadj"),
              pybind11::arg("adjcwgt"),
              pybind11::arg("adjncy"),
              pybind11::arg("suppress_output") = true,
              pybind11::arg("seed") = 0,
              pybind11::arg("time_limit") = 1.0,
              pybind11::arg("cluster_upperbound") = 0,
              pybind11::arg("leiden") = false,
              pybind11::arg("leiden_theta") = 0.01);
}
