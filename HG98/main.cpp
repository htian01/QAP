#include <stdlib.h>
#include <stdio.h>
#include <omp.h> 

#include "dataset.h"
#include "DP.h"

using namespace std;


int main(int argc, char* argv[])
{
    
    string load_path = argv[1];
    string save_path = argv[2];
    string problem_name = argv[3];
    int max_steps = stoi(argv[4]);

    string load_file_name = load_path + "/" + problem_name + ".dat";
    string save_file_name = save_path + "/" + problem_name + "_reparam.dd";

    float* C;
    int N;
    //read .dat problem
    tie(C, N) = read_QAPLIB(load_file_name);
    printf("%s\n", problem_name.c_str());
    if (N < omp_get_max_threads())
    {
        omp_set_num_threads(N);
    }

    auto dp = DP(C, N, max_steps);
    auto bound = dp.solve();
    auto cost = dp.get_cost(dp.results);
    free(C);

    // get_reparametrized_matrix 
    float* reparametrized_C = (float*)calloc(N * N * N * N, sizeof(float));
    dp.get_reparametrized_matrix(reparametrized_C);
    Model* reparametrized_m = matrix_to_model(reparametrized_C, N);
    reparametrized_m->save_dd_format(save_file_name.c_str(), bound, cost);
}
