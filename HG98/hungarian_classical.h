#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <stack>
#include <vector>

class Hungarian_Classical
{
public:
	Hungarian_Classical(float* matrix_2d, const int size);
	~Hungarian_Classical();
	float solve();
	float* C;
	float* ori_C;
	const int N;
	int* Ar;
	float cost;
private:
	int* Ac;
	int* Vr;
	int* Vc;
	int* Pr;
	int* Pc;
	int match_count;
	void set_var_sentinels();
	bool check();
	bool search();
	void augment(int j);
	void update();

};
float min_with_stride(float* array, const int num_elements, const int stride = 1)
{
	float min_val = array[0];
	for (int i = 1; i < num_elements; i++)
	{
		min_val = std::min(min_val, array[i * stride]);
	}
	return min_val;
}

Hungarian_Classical::Hungarian_Classical(float* matrix_2d, const int size) : C(matrix_2d), N(size)
{
	cost = 0;
	match_count = 0;
	ori_C = (float*)malloc(N * N * sizeof(float));
	for (int i = 0; i < N * N; i++)
	{
		ori_C[i] = C[i];
	}

	Ar = (int*)malloc(N * sizeof(int));
	Ac = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)
	{
		Ar[i] = -1;
		Ac[i] = -1;
	}
	Vr = (int*)malloc(N * sizeof(int));
	Vc = (int*)malloc(N * sizeof(int));
	Pr = (int*)malloc(N * sizeof(int));
	Pc = (int*)malloc(N * sizeof(int));
	set_var_sentinels();
}

Hungarian_Classical::~Hungarian_Classical()
{
	free(ori_C);
	free(Ar);
	free(Ac);
	free(Vr);
	free(Vc);
	free(Pr);
	free(Pc);
}

float Hungarian_Classical::solve()
{
	// row reduce
	for (int i = 0; i < N; i++)
	{
		float Dr_i = min_with_stride(C + i * N, N, 1);
		for (int j = 0; j < N; j++)
		{
			C[i * N + j] -= Dr_i;
		}
	}
	// col reduce
	for (int j = 0; j < N; j++)
	{
		float Dc_j = min_with_stride(C + j, N, N);
		for (int i = 0; i < N; i++)
		{
			C[i * N + j] -= Dc_j;
		}
	}
	check();
	while (true)
	{
		if (!search())
		{
			if (check())
			{
				for (int i = 0; i < N; i++)
				{
					int j = Ar[i];
					cost += ori_C[i * N + j];
				}
				return cost;
			}
		}
	}
}

void Hungarian_Classical::set_var_sentinels()
{
	for (int i = 0; i < N; i++)
	{
		Vr[i] = 0;
		Vc[i] = 0;
		Pr[i] = -1;
		Pc[i] = -1;
	}
}

bool Hungarian_Classical::check()
{
	match_count = 0;
	set_var_sentinels();
	for (int i = 0; i < N; i++)
	{
		if (Ar[i] != -1)
		{
			Vr[i] = 1;
			match_count += 1;
		}
	}
	if (match_count == N)
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool Hungarian_Classical::search()
{
	std::stack<int> slack;
	std::vector<std::stack<int>> Z(N);
	for (int i = 0; i < N; i++)
	{
		if (Vr[i] == 0)
		{
			slack.push(i);
			for (int j = 0; j < N; j++)
			{
				if (C[i * N + j] == 0)
				{
					Z[i].push(j);
				}
			}
		}
	}
	while (!slack.empty())
	{
		int i = slack.top();
		slack.pop();
		while (!Z[i].empty())
		{
			int j = Z[i].top();
			Z[i].pop();
			int i_new = Ac[j];
			if (i_new == i)
			{
				continue;
			}
			if (Vc[j] == 0)
			{
				Pc[j] = i;
				if (i_new == -1)
				{
					augment(j);
					return false;
				}
				else
				{
					slack.push(i_new);
					Pr[i_new] = j;
					Vr[i_new] = 0;
					Vc[j] = 1;
				}
			}
		}
	}
	update();
	return true;
}


void Hungarian_Classical::augment(int j)
{
	int c_cur = j;
	int r_cur = -1;
	while (c_cur != -1)
	{
		r_cur = Pc[c_cur];
		Ar[r_cur] = c_cur;
		Ac[c_cur] = r_cur;
		c_cur = Pr[r_cur];
	}
}

void Hungarian_Classical::update()
{
	float theta = 1e10;
	for (int i = 0; i < N; i++)
	{
		if (Vr[i] != 0)
		{
			continue;
		}
		for (int j = 0; j < N; j++)
		{
			if (Vc[j] != 0)
			{
				continue;
			}
			else
			{
				theta = std::min(theta, C[i * N + j]);
			}
		}
	}

	for (int k = 0; k < N; k++)
	{
		if (Vr[k] == 0)
		{
			for (int n = 0; n < N; n++)
			{
				C[k * N + n] -= theta;
			}
		}
		if (Vc[k] == 1)
		{
			for (int n = 0; n < N; n++)
			{
				C[n * N + k] += theta;
			}
		}
	}
}
