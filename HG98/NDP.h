#pragma once
#include "ODP.h"
class NDP :
	public ODP
{
	using ODP::ODP;
	void move_complementary(const int i, const int j, const int row, const int col) override;
};
void NDP::move_complementary(const int i, const int j, const int row, const int col)
{
	int k = row + int(row >= i);
	int n = col + int(col >= j);
	int new_row = i - int(i >= k);
	int new_col = j - int(j >= n);
	// self.sub_matrix[k][n][new_row][new_col] += self.sub_matrix[i][j][row][col]
	// self.sub_matrix[i][j][row][col] = 0

	float complementary_sum = sub_matrix[SUBMATRIX_INDEX(k, n, new_row, new_col, N)] + sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)];
	sub_matrix[SUBMATRIX_INDEX(k, n, new_row, new_col, N)] = complementary_sum / 2.;
	sub_matrix[SUBMATRIX_INDEX(i, j, row, col, N)] = complementary_sum / 2.;
}