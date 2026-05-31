/******************************************************************************
 * parse_parameters.h
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/


#ifndef PARSE_PARAMETERS_GPJMGSM8
#define PARSE_PARAMETERS_GPJMGSM8

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_MPI
#include <mpi.h>
#else
#include "tools/pseudo_mpi.h"
#endif
#include <sstream>
#include "configuration.h"
#include "clustering/leiden_config.h"

static LeidenConfig _default_leiden_config;

int parse_parameters(int argn, char **argv,
                     PartitionConfig & partition_config,
                     std::string & graph_filename,
                     bool & is_graph_weighted,
                     bool & suppress_program_output,
                     bool & recursive,
                     LeidenConfig & leiden_config = _default_leiden_config) {

        const char *progname = argv[0];

        // Setup argtable parameters.
        struct arg_lit *help                                 = arg_lit0(NULL, "help","Print help.");
        struct arg_dbl *mh_mutate_fraction                   = arg_dbl0(NULL, "mh_mutate_fraction", NULL, "Fraction to determine the number of clusters that will be split by the mutation operator.");
        struct arg_int *local_partitioning_repetitions       = arg_int0(NULL, "local_partitioning_repetitions", NULL, "Number of local repetitions.");
        struct arg_str *filename                             = arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file to partition.");
        struct arg_str *filename_output                      = arg_str0(NULL, "output_filename", NULL, "Specify the name of the output file (that contains the partition).");
        struct arg_int *user_seed                            = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
        struct arg_str *input_partition                      = arg_str0(NULL, "input_partition", NULL, "Input partition to use.");
        struct arg_dbl *time_limit                           = arg_dbl1(NULL, "time_limit", NULL, "Time limit in s. Default 0s .");
        struct arg_lit *mh_print_log                         = arg_lit0(NULL, "mh_print_log", "Each PE prints a logfile (timestamp, edgecut).");
        struct arg_int *cluster_upperbound                   = arg_int0(NULL, "cluster_upperbound", NULL, "Set a size-constraint on the size of a cluster. Default: none");
        struct arg_int *label_propagation_iterations         = arg_int0(NULL, "label_propagation_iterations", NULL, "Set the number of label propgation iterations. Default: 10.");

        // for graph clustering we need some own parameters (by BSc)
        struct arg_dbl *lm_minimum_quality_improvement              = arg_dbl0(NULL, "lm_minimum_quality_improvement", NULL, "Minimum quality (modularity) improvement necessary to perform another turn in the Louvain method. Default: 1e-6.");
        struct arg_int *lm_number_of_label_propagation_levels       = arg_int0(NULL, "lm_number_of_label_propagation_levels", NULL, "Number of label propagation levels used before the Louvain method. Default: 0.");
        struct arg_int *lm_number_of_label_propagation_iterations   = arg_int0(NULL, "lm_number_of_label_propagation_iterations", NULL, "Number of iterations per label propagation level. Default: 3.");
        struct arg_int *lm_cluster_coarsening_factor                = arg_int0(NULL, "lm_cluster_coarsening_factor", NULL, "Factor relative to the number of nodes that limits the maximum cluster size for size constrained label propagation. If this factor is 0 or 1, then no size constraint is used. Default: 0.");
        struct arg_str *output_log_json                             = arg_str0(NULL, "output_log_json", NULL, "File to write the log in JSON format into it. Use \"-\" to write to STDOUT. (Default: disabled)");

        struct arg_int *mh_pool_size                         = arg_int0(NULL, "mh_pool_size", NULL, "MetaHeuristic Pool Size.");
        struct arg_lit *leiden                               = arg_lit0(NULL, "leiden", "Use Leiden algorithm (guarantees connected communities).");
        struct arg_dbl *leiden_theta                         = arg_dbl0(NULL, "leiden_theta", NULL, "Leiden refinement temperature parameter. Default: 0.01.");
        struct arg_end *end                                  = arg_end(100);

        // Define argtable.
        void* argtable[] = {
                help, filename, user_seed,
#ifdef MODE_KAFFPAE
                time_limit,
                //mh_pool_size,
                //mh_mutate_fraction,
		//mh_print_log,
                //local_partitioning_repetitions,
                //input_partition,
                filename_output,
                leiden,
                leiden_theta,
#elif defined MODE_EVALUATOR
                input_partition,
#elif defined MODE_CLUSTERING
    // for graph clustering we need some own parameters (by BSc)
    lm_minimum_quality_improvement,
    lm_number_of_label_propagation_levels,
    lm_number_of_label_propagation_iterations,
    // if we use size constrained label propagation
    lm_cluster_coarsening_factor,
    output_log_json,
    leiden,
    leiden_theta,
#endif
                end
        };
        // Parse arguments.
        int nerrors = arg_parse(argn, argv, argtable);

        // Catch case that help was requested.
        if (help->count > 0) {
                int rank, size;
                MPI_Comm communicator = MPI_COMM_WORLD; 
                MPI_Comm_rank( communicator, &rank);
                MPI_Comm_size( communicator, &size);

                if( rank == 0 ) {
                printf("Usage: %s", progname);
                arg_print_syntax(stdout, argtable, "\n");
                arg_print_glossary(stdout, argtable,"  %-40s %s\n");
                arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
                }
                MPI_Finalize();
                return 1;
        }


        if (nerrors > 0) {
                arg_print_errors(stderr, end, progname);
                printf("Try '%s --help' for more information.\n",progname);
                arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
                return 1;
        }

        if(filename->count > 0) {
                graph_filename = filename->sval[0];
        }

        partition_config.k = 2;
        recursive = false;

        configuration cfg;
        cfg.standard(partition_config);

        cfg.strong(partition_config);

        
        if(filename_output->count > 0) {
                partition_config.filename_output = filename_output->sval[0];
        }

        if(mh_print_log->count > 0) {
                partition_config.mh_print_log = true;
        }

        if(time_limit->count > 0) {
                partition_config.time_limit = time_limit->dval[0];
        }

        if(mh_pool_size->count > 0) {
                partition_config.mh_pool_size = mh_pool_size->ival[0];
        }

        if(input_partition->count > 0) {
                partition_config.input_partition = input_partition->sval[0];
        }

        if (label_propagation_iterations->count > 0) {
                partition_config.label_iterations = label_propagation_iterations->ival[0];
        }

        if (cluster_upperbound->count > 0) {
                partition_config.cluster_upperbound = cluster_upperbound->ival[0];
        } else {
                partition_config.cluster_upperbound = std::numeric_limits< NodeWeight >::max()/2;
        }

        // for graph clustering we need some own parameters (by BSc)
        if (lm_minimum_quality_improvement->count > 0) {
            partition_config.lm_minimum_quality_improvement = lm_minimum_quality_improvement->dval[0];
        }

        if (lm_number_of_label_propagation_levels->count > 0) {
            partition_config.lm_number_of_label_propagation_levels = static_cast<unsigned>(lm_number_of_label_propagation_levels->ival[0]);
        }

        if (lm_number_of_label_propagation_iterations->count > 0) {
            partition_config.lm_number_of_label_propagation_iterations = static_cast<unsigned>(lm_number_of_label_propagation_iterations->ival[0]);
        }

        if (lm_cluster_coarsening_factor->count > 0) {
            partition_config.lm_cluster_coarsening_factor = static_cast<unsigned>(lm_cluster_coarsening_factor->ival[0]);
        }

        if (output_log_json->count > 0) {
            partition_config.outputLogJsonFileName = output_log_json->sval[0];
        }

        if(leiden->count > 0) {
                leiden_config.enabled = true;
        }

        if(leiden_theta->count > 0) {
                leiden_config.theta = leiden_theta->dval[0];
        }

        if(mh_mutate_fraction->count > 0) {
                partition_config.mh_mutate_fraction = mh_mutate_fraction->dval[0];
        }

        if(local_partitioning_repetitions->count > 0) {
                partition_config.local_partitioning_repetitions = local_partitioning_repetitions->ival[0];
        }


        return 0;
}

#endif /* end of include guard: PARSE_PARAMETERS_GPJMGSM8 */
