namespace legalizer
{

enum class LGSortingPolicy { RANDOM, SIZE, OVERFLOW };
enum class LGSortingAlgorithm { S_MERGE, M_MERGE, S_BITONIC, M_BITONIC, CUB };

class LGHyperParameters
{
  public:
    bool flag_iccad17;
    bool standard_admm;
    int log_freq;
    int num_cpu_threads;
    int max_disp_in_row_height;
    int max_admm_iter;
    int admm_terminate_thr;
    int x_hint;
    int y_hint;
    int iter_update_partition;
    float max_disp_coeff;
    float tech_penalty;
    float edge_spacing_penalty;

    // LEGALM style ADMM
    int iter_threshold_dual_update;
    int iter_update_dual_step;
    int iter_perturbation_on;
    float coeff_to_init_rho;
    float coeff_to_init_dual;
    float coeff_to_update_dual_step;

    // QP ADMM
    bool qp_admm_on;
    bool qp_admm_dp_site_align;
    int qp_admm_max_iter;
    int qp_admm_min_iter;
    float qp_admm_rho;
    
    // Adaptive Search
    bool use_adaptive_search;
    int stagnation_window;
    int min_iter_to_call_adaptive_search;
    int cooltime_to_call_adaptive_search;
    float stagnation_threshold;
    float ovf_ratio_to_call_adaptive_search;
    float coeff_to_neglect_displace;

    // Sorting Policy
    int random_sorting_seed;
    LGSortingPolicy sorting_policy;
    LGSortingAlgorithm sorting_algorithm;
};

}
