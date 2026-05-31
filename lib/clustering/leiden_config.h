/******************************************************************************
 * leiden_config.h
 *
 * Source of VieClus -- Vienna Graph Clustering
 *****************************************************************************/


#ifndef LEIDEN_CONFIG_H
#define LEIDEN_CONFIG_H

struct LeidenConfig {
        bool enabled;
        double theta;
        LeidenConfig() : enabled(false), theta(0.01) {}
        LeidenConfig(bool e, double t) : enabled(e), theta(t) {}
};

#endif // LEIDEN_CONFIG_H
