#pragma once
#include "NDP.h"
class DP:
	public NDP
{
public:
	using NDP::NDP;
	float tmp;
protected:
	void step3() override;
};

void DP::step3()
{
	float return_super_leader = super_leader * exp(-float(steps) / (float(max_steps) / 20.0));
	super_leader -= return_super_leader;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			leader_matrix[i * N + j] += return_super_leader / N;
		}
	}
	NDP::step3();
}