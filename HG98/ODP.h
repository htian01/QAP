#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <chrono> 
#include <omp.h> 
#define PRECISION 1e-6
#include "hungarian_classical.h"

using namespace std;
using namespace chrono;
using clock_type = std::chrono::high_resolution_clock;

void print_2d_matrix(float* x, const int N);
inline int SUBMATRIX_INDEX(int i, int j, int k, int n, int N);
class ODP
{
public:
	ODP(float* matrix_4d, const int size, const int max_steps = 100);
	virtual ~ODP();
	virtual float solve();
	void get_reparametrized_matrix(float* new_C);
	float get_cost(int* results);
	float* C;  // [N,N,N,N]
	const int N;
	float super_leader;
	int* results; // [N], results[i] = j means x_{ij} = 1;
	float* sub_matrix;  // [N,N,N-1,N-1]
	float* leader_matrix;  //[N,N]
	float* step2_leader_matrix;  //[N,N]
	int steps;
	const int max_steps;
protected:
	static void build_sub_matrix(float* cost_matrix, const int N, float* sub_matrix);
	static void build_leader_matrix(float* cost_matrix, const int N, float* leader_matrix);
	void upper_most_move();
	void solve_sub_matrix(const int i, const int j);
	bool check_solution(int* results);
	void move_sub_matrix_complementary_to_lower(const int i, const int j);

	virtual void move_complementary(const int i, const int j, const int row, const int col);
	virtual void step1();
	virtual bool step2();
	virtual void step3();
	virtual void step4();
};

inline int SUBMATRIX_INDEX(int i, int j, int k, int n, int N)
{
	return (i * (N * (N - 1) * (N - 1)) + j * ((N - 1) * (N - 1)) + k * (N - 1) + n);
}
//#define SUBMATRIX_INDEX(i,j,k,n,N) (i*(N*(N-1)*(N-1)) + j * ((N-1)*(N-1)) + k*(N-1) + n)

void print_2d_matrix(float* x, const int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%.10f,", x[i * N + j]);
		}
		printf("\n");
	}
}


ODP::ODP(float* matrix_4d, const int size, const int max_steps) : C(matrix_4d), N(size), max_steps(max_steps)
{
	steps = 0;
	super_leader = 0;

	sub_matrix = (float*)malloc(N * N * (N - 1) * (N - 1) * sizeof(float));
	build_sub_matrix(C, N, sub_matrix);
	leader_matrix = (float*)malloc(N * N * sizeof(float));
	build_leader_matrix(C, N, leader_matrix);
	results = (int*)malloc(N * sizeof(int));
	step2_leader_matrix = (float*)malloc(N * N * sizeof(float));
}


ODP::~ODP()
{
	free(results);
	free(sub_matrix);
	free(leader_matrix);
	free(step2_leader_matrix);
}

float ODP::solve()
{
	auto start = clock_type::now();
	step1();
	while (true)
	{
		bool is_continue = step2();
		auto during = std::chrono::duration_cast<std::chrono::duration<double>>(clock_type::now() - start).count();
		printf("it=%d\t lb=%.4lf\t t=%.4lf\n", steps - 1, super_leader, during);
		if (is_continue)
		{
			step3();
			step4();
		}
		else
		{
			return super_leader;
		}
	}
	return super_leader;
}

void ODP::get_reparametrized_matrix(float* new_C)
{
	step3();
	upper_most_move();
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int row = 0; row < N - 1; row++)
			{
				for (int col = 0; col < N - 1; col++)
				{
					//scost_matrix[i][j][row + (row >= i)][col + (col >= j)] = sub_matrix[i][j][row][col] 
					new_C[i * N * N * N + j * N * N + (row + int(row >= i)) * N + (col + int(col >= j))] = sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)];
				}
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// cost_matrix[i][j][i][j] = leader_matrix[i][j]
			new_C[i * N * N * N + j * N * N + i * N + j] = leader_matrix[i * N + j];
		}
	}
}

float ODP::get_cost(int* results)
{
	float cost = 0;
	for (int i = 0; i < N; i++)
	{
		int j = results[i];
		for (int k = 0; k < N; k++)
		{
			int n = results[k];
			cost += C[i * N * N * N + j * N * N + k * N + n];
		}
	}
	return cost;
}

void ODP::build_sub_matrix(float* cost_matrix, const int N, float* sub_matrix)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int row = 0; row < N - 1; row++)
			{
				for (int col = 0; col < N - 1; col++)
				{
					// sub_matrix[i][j][row][col] = cost_matrix[i][j][row + (row >= i)][col + (col >= j)]
					sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)] = cost_matrix[i * N * N * N + j * N * N + (row + int(row >= i)) * N + (col + int(col >= j))];
				}
			}
		}
	}
}

void ODP::build_leader_matrix(float* cost_matrix, const int N, float* leader_matrix)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// leader_matrix[i][j] = cost_matrix[i][j][i][j]
			leader_matrix[i * N + j] = cost_matrix[i * N * N * N + j * N * N + i * N + j];
		}
	}
}

void ODP::move_complementary(const int i, const int j, const int row, const int col)
{
	int k = row + int(row >= i);
	int n = col + int(col >= j);
	int new_row = i - int(i >= k);
	int new_col = j - int(j >= n);
	// self.sub_matrix[k][n][new_row][new_col] += self.sub_matrix[i][j][row][col]
	// self.sub_matrix[i][j][row][col] = 0
	sub_matrix[SUBMATRIX_INDEX(k, n, new_row, new_col, N)] += sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)];
	sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)] = 0;

}

void ODP::upper_most_move()
{
#pragma omp parallel for
	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int row = 0; row < i; row++)
			{
				for (int col = 0; col < N - 1; col++)
				{
					move_complementary(i, j, row, col);
				}
			}
		}
	}
}

void ODP::solve_sub_matrix(const int i, const int j)
{
	auto cost = Hungarian_Classical(sub_matrix + SUBMATRIX_INDEX(i, j, 0, 0, N), N - 1).solve();
	leader_matrix[i * N + j] += cost;
}

bool ODP::check_solution(int* results)
{
	for (int i = 0; i < N; i++)
	{
		int j = results[i];
		for (int k = 0; k < N; k++)
		{
			int n = results[k];
			if (i != k && j != n)
			{
				int row = k - int(k >= i);
				int col = n - int(n >= j);
				if (sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)] > PRECISION)
					//if (sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)] != 0)
				{
					return false;
				}
			}
		}
	}
	return true;
}

void ODP::move_sub_matrix_complementary_to_lower(const int i, const int j)
{
	for (int row = i; row < N - 1; row++)
	{
		for (int col = 0; col < N - 1; col++)
		{
			move_complementary(i, j, row, col);
		}
	}
}

void ODP::step1()
{
	upper_most_move();
	for (int i = 0; i < N; i++)
	{
#pragma omp parallel for
		for (int j = 0; j < N; j++)
		{
			solve_sub_matrix(i, j);
			move_sub_matrix_complementary_to_lower(i, j);
		}
	}
	upper_most_move();
}

bool ODP::step2()
{
	steps += 1;
	auto hungarian = Hungarian_Classical(leader_matrix, N);
	auto cost = hungarian.solve();
	memcpy(step2_leader_matrix, leader_matrix, N * N * sizeof(float));
	super_leader += cost;
	//printf("Step %d Bound %.3f Cost %.3f \n", steps, super_leader, get_cost(hungarian.Ar));
	if (check_solution(hungarian.Ar))
	{
		memcpy(results, hungarian.Ar, N * sizeof(int));
		return false;
	}
	else if (steps >= max_steps)
	{
		memcpy(results, hungarian.Ar, N * sizeof(int));
		return false;
	}
	else if (cost > PRECISION)
		//else if (cost > 0)
	{
		return true;
	}
	else
	{
		memcpy(results, hungarian.Ar, N * sizeof(int));
		return false;
	}
}

void ODP::step3()
{
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			float num = leader_matrix[i * N + j];
			leader_matrix[i * N + j] = 0;
			const int div = N - 1;
			float distribution_number = num / float(div);
			for (int row = 0; row < N - 1; row++)
			{
				for (int col = 0; col < N - 1; col++)
				{
					sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)] += distribution_number;
				}
			}
		}
	}
	upper_most_move();
}

void ODP::step4()
{
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (step2_leader_matrix[i * N + j] == 0)
			{
				solve_sub_matrix(i, j);
			}
		}

	}
	upper_most_move();  // comment this lead to little worese bound but boost the speed
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (step2_leader_matrix[i * N + j] != 0)
			{
				solve_sub_matrix(i, j);
			}
		}
	}
}
