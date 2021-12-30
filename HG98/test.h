#include <stdlib.h>
#include <stdio.h>
#include <windows.h>

#include "dataset.h"
#include "hungarian_classical.h"
#include "DP.h"

using namespace std;

void test_hungarian()
{
    const int N = 12;
    float cost_matrix[] = { 0.36, 0.33, 0.00, 0.38, 0.07, 0.45, 0.63, 0.67, 1.17, 0.64, 3.42, 0.54,
0.61, 0.74, 0.94, 1.05, 1.05, 0.92, 0.19, 1.13, 0.86, 1.03, 0.00, 1.47,
1.43, 1.43, 1.54, 1.42, 1.75, 1.31, 1.74, 1.51, 1.37, 0.00, 1.35, 1.70,
1.89, 1.85, 2.20, 2.36, 1.98, 2.02, 2.15, 2.15, 3.18, 2.07, 6.05, 0.00,
3.05, 3.15, 2.95, 2.99, 2.80, 0.00, 3.12, 3.09, 3.81, 2.98, 3.49, 3.28,
3.51, 4.51, 3.68, 3.54, 5.88, 3.78, 0.00, 3.83, 3.87, 3.57, 3.57, 5.30,
4.30, 4.57, 4.82, 4.39, 6.74, 4.59, 4.38, 0.00, 4.38, 4.44, 4.17, 5.54,
0.00, 4.76, 4.85, 5.01, 5.02, 5.03, 4.91, 4.96, 0.00, 4.69, 5.50, 5.20,
5.41, 5.60, 0.26, 0.00, 5.35, 5.73, 5.61, 5.61, 5.75, 5.65, 7.13, 0.00,
5.98, 6.35, 6.13, 6.21, 6.03, 6.18, 6.26, 6.47, 0.00, 6.29, 6.25, 6.48,
6.81, 0.00, 6.77, 6.83, 6.66, 7.18, 6.74, 6.83, 6.77, 6.89, 6.72, 7.01,
7.50, 7.57, 7.78, 0.00, 0.00, 7.57, 7.47, 7.68, 7.73, 7.66, 9.93, 7.80, };
    auto hungarian = Hungarian_Classical(cost_matrix, N);
    auto cost = hungarian.solve();
    printf("Solution:"); // (0, 2), (1, 1), (2, 0)
    for (int i = 0; i < N; i++)
    {
        printf(" (%d, %d) ", i, hungarian.Ar[i]);
    }
    printf("Cost: %.2f", cost); // 407
}

void test_DP()
{
    float C[] = {
        105,0,0,0, 0,60,120,30, 0,70,140,35, 0,20,40,10,
        0, 105, 0, 0, 108, 0, 54, 24, 126, 0, 63, 28, 36, 0, 18, 8,
        0, 0, 105, 0, 30, 36, 0, 48, 35, 42, 0, 56, 10, 12, 0, 16,
        0, 0, 0, 105, 48, 0, 90, 0, 56, 0, 105, 0, 16, 0, 30, 0,
        0, 60, 120, 30, 90, 0, 0, 0, 0, 50, 100, 25, 0, 60, 120, 30,
        108, 0, 54, 24, 0, 90, 0, 0, 90, 0, 45, 20, 108, 0, 54, 24,
        30, 36, 0, 48, 0, 0, 90, 0, 25, 30, 0, 40, 30, 36, 0, 48,
        48, 0, 90, 0, 0, 0, 0, 90, 40, 0, 75, 0, 48, 0, 90, 0,
        0, 70, 140, 35, 0, 50, 100, 25, 105, 0, 0, 0, 0, 10, 20, 5,
        126, 0, 63, 28, 90, 0, 45, 20, 0, 105, 0, 0, 18, 0, 9, 4,
        35, 42, 0, 56, 25, 30, 0, 40, 0, 0, 105, 0, 5, 6, 0, 8,
        56, 0, 105, 0, 40, 0, 75, 0, 0, 0, 0, 105, 8, 0, 15, 0,
        0, 20, 40, 10, 0, 60, 120, 30, 0, 10, 20, 5, 90, 0, 0, 0,
        36, 0, 18, 8, 108, 0, 54, 24, 18, 0, 9, 4, 0, 90, 0, 0,
        10, 12, 0, 16, 30, 36, 0, 48, 5, 6, 0, 8, 0, 0, 90, 0,
        16, 0, 30, 0, 48, 0, 90, 0, 8, 0, 15, 0, 0, 0, 0, 90,
    };
    const int N = 4;
    LARGE_INTEGER t1, t2, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    auto dp = DP(C, 4, 2000);
    auto cost = dp.solve();
    QueryPerformanceCounter(&t2);
    auto time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart * 1000;
    printf("Time takes: %lf ms\n", time);
    printf("Cost: %.2f", cost); // 793
}

void test_DP_2()
{
    const int N = 4;
    float D[] = {
        0, 22, 53,53,
        22,0,40,62,
        53,40,0,55,
        53,62,55,0,
    };
    float F[] = {
        0,3,0,2,
        3,0,0,1,
        0,0,0,4,
        2,1,4,0,
    };

    float* C = (float*)malloc(N * N * N * N * sizeof(float));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                for (int n = 0; n < N; n++)
                {
                    C[i * N * N * N + j * N * N + k * N + n] = D[j * N + n] * F[i * N + k];
                }
            }
        }
    }

    LARGE_INTEGER t1, t2, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    auto dp = DP(C, 4, 2000);
    auto cost = dp.solve();
    QueryPerformanceCounter(&t2);
    auto time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart * 1000;
    printf("Time takes: %lf ms\n", time);
    printf("Cost: %.2f", cost); // 793
}


void test_QAPLIB(const string file_name, const int max_steps=100)
{
    float* C;
    int N;
    tie(C, N) = read_QAPLIB(file_name);

    LARGE_INTEGER t1, t2, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    auto dp = DP(C, N, max_steps);
    auto bound = dp.solve();
    auto cost = dp.get_cost(dp.results);
    QueryPerformanceCounter(&t2);
    auto time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    printf("For task %s N: %d Steps: %d Time takes: %.2lf s Lower Bound: %.2f Cost: %.2f\n", file_name.c_str(), N, dp.steps, time, bound, cost);
    free(C);
    printf("Solution: (");
    for (int i = 0; i < N; i++)
    {
        printf("%d, ", dp.results[i]);
    }
    printf(")\n");


}

void test_dd_format_matrix_conversion()
{
    string file_name = "D:\\dataset\\QAP\\datasets\\datasets\\car-motor\\car1.dd";
    const auto max_steps = 100;
    float* C;
    int N;
    tie(C, N) = read_dd_format(file_name);
    Model* new_m = matrix_to_model(C, N);
    float* new_C;
    tie(new_C, N) = model_to_matrix(new_m);

    const char* save_file_name = "D:\\dataset\\QAP\\datasets\\datasets\\car-motor\\car1_saved.dd";
    new_m->save_dd_format(save_file_name);
    float* saved_C;
    tie(saved_C, N) = read_dd_format(file_name);

    float diff = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            for (size_t k = 0; k < N; k++)
            {
                for (size_t p = 0; p < N; p++)
                {
                    if (p != k)
                    {
                        diff += C[i * N * N * N + j * N * N + k * N + p] - new_C[i * N * N * N + j * N * N + k * N + p];
                        diff += C[k * N * N * N + p * N * N + i * N + j] - new_C[k * N * N * N + p * N * N + i * N + j];
                    }
                }
            }
        }
    }
    printf("Diff: %.3f\n", diff);
    diff = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            for (size_t k = 0; k < N; k++)
            {
                for (size_t p = 0; p < N; p++)
                {
                    if (p != k)
                    {
                        diff += C[i * N * N * N + j * N * N + k * N + p] - saved_C[i * N * N * N + j * N * N + k * N + p];
                        diff += C[k * N * N * N + p * N * N + i * N + j] - saved_C[k * N * N * N + p * N * N + i * N + j];
                    }
                }
            }
        }
    }
    printf("Diff: %.3f\n", diff);

    for (int i = 0; i < N * N * N * N; i++)
    {
        C[i] += 1;
    }
    for (int i = 0; i < N * N * N * N; i++)
    {
        new_C[i] += 1;
    }
    for (int i = 0; i < N * N * N * N; i++)
    {
        saved_C[i] += 1;
    }
    LARGE_INTEGER t1, t2, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    auto dp = DP(C, N, max_steps);
    auto bound = dp.solve() - N * N;
    auto cost = dp.get_cost(dp.results) - N * N;
    QueryPerformanceCounter(&t2);
    auto time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    free(C);
    printf("For task %s N: %d Steps: %d Time takes: %.2lf s Lower Bound: %.2f Cost: %.2f\n", file_name.c_str(), N, dp.steps, time, bound, cost);
    printf("Solution: (");
    for (int i = 0; i < N; i++)
    {
        printf("%d, ", dp.results[i]);
    }
    printf(")\n");
}
void bench_QAPLIB()
{
    test_QAPLIB("./qaplib/had16.dat", 2000);
    /*
    test_QAPLIB("./qaplib/had18.dat", 2000);
    test_QAPLIB("./qaplib/had20.dat", 2000);
    test_QAPLIB("./qaplib/kra30a.dat", 2000);
    test_QAPLIB("./qaplib/kra30b.dat", 2000);
    test_QAPLIB("./qaplib/nug12.dat", 2000);
    test_QAPLIB("./qaplib/nug15.dat", 2000);
    test_QAPLIB("./qaplib/nug30.dat", 2000);
    */
}

void test_car_motor(const string file_name, const int max_steps = 100)
{
    float* C;
    int N;
    tie(C, N) = read_dd_format(file_name);
    for (int i = 0; i < N * N * N * N; i++)
    {
        C[i] += 1;
    }

    LARGE_INTEGER t1, t2, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    auto dp = DP(C, N, max_steps);
    auto bound = dp.solve() - N * N;
    auto cost = dp.get_cost(dp.results) - N * N;
    QueryPerformanceCounter(&t2);
    auto time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    free(C);
    printf("For task %s N: %d Steps: %d Time takes: %.2lf s Lower Bound: %.2f Cost: %.2f\n", file_name.c_str(), N, dp.steps, time, bound, cost);
    printf("Solution: (");
    for (int i = 0; i < N; i++)
    {
        printf("%d, ", dp.results[i]);
    }
    printf(")\n");
}

void bench_car_motor()
{
    for (int i = 1; i <= 30; i++)
    {
        string file_name = string("../../dataset/car_motor/car/car") + to_string(i) + string(".txt");
        test_car_motor(file_name, 500);
    }
    for (int i = 1; i <= 20; i++)
    {
        string file_name = string("../../dataset/car_motor/motor/motor") + to_string(i) + string(".txt");
        test_car_motor(file_name, 500);
    }
}

void reparametrized_QAPLIB(const string problem_name)
{
    const string load_path = "./qaplib/";
    const string save_path = "./qaplib/";

    const string load_file_name = load_path + problem_name + ".dat";
    const string save_file_name = save_path + problem_name + "_reparametrized.dd";
    const int max_steps = 2000;

    float* C;
    int N;
    tie(C, N) = read_QAPLIB(load_path + problem_name + ".dat");
    LARGE_INTEGER t1, t2, tc;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    auto dp = DP(C, N, max_steps);
    auto bound = dp.solve();
    auto cost = dp.get_cost(dp.results);
    QueryPerformanceCounter(&t2);
    auto time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    printf("For task %s N: %d Steps: %d Time takes: %.2lf s Lower Bound: %.2f Cost: %.2f\n", problem_name.c_str(), N, dp.steps, time, bound, cost);
    free(C);
    printf("Solution: (");
    for (int i = 0; i < N; i++)
    {
        printf("%d, ", dp.results[i]);
    }
    printf(")\n");

    // get_reparametrized_matrix
    float* reparametrized_C = (float*)calloc(N * N * N * N, sizeof(float));
    dp.get_reparametrized_matrix(reparametrized_C);
    Model* reparametrized_m = matrix_to_model(reparametrized_C, N);
    float max_cost = 0;
    for (int i = 0; i < N * N * N * N; i++)
    {
        max_cost = max(max_cost, reparametrized_C[i]);
    }
    //int incentive = (int(N * max_cost) / 100 + 1) * 100;
    //printf("Incentive %d", incentive);
    /*
    for (auto it = reparametrized_m->assignments.begin(); it != reparametrized_m->assignments.end(); it++)
    {
        assert(get<2>(*it) == 0);
        get<2>(*it) -= incentive;
    }
    */
    reparametrized_m->save_dd_format(save_file_name.c_str(), bound, cost);
}
