// Copyright (c) 2015-2017, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause
#include "assert.h"
#include "src/core/matrix.hpp"

using namespace raptor;

void add(CSRMatrix* B, CSRMatrix* C, double alpha, bool duplicates)
{
    int start, end;
    int C_nnz = nnz + B->nnz;
    C->cols.resize(C_nnz);
    C->data.resize(C_nnz);

    C_nnz = 0;
    C->rowptr[0] = 0;
    for (int i = 0; i < n_rows; i++)
    {
        start = rowptr[i];
        end = rowptr[i+1];
        std::copy(cols.begin() + start,
                cols.begin() + end,
                C->cols.begin() + C_nnz);
        if (data.size())
        {
            std::copy(data.begin() + start,
                    data.begin() + end,
                    C->data.begin() + C_nnz);
        }
        C_nnz += (end - start);

        start = B->rowptr[i];
        end = B->rowptr[i+1];
        std::copy(B->cols.begin() + start,
                B->cols.begin() + end,
                C->cols.begin() + C_nnz);
        if (B->data.size())
        {
            for (int j = 0; j < (end - start); j++)
            {
                C->data[C_nnz+j] = alpha*B->data[start+j];
            }
        }
        C_nnz += (end - start);

        C->rowptr[i+1] = C_nnz;
    }
    C->nnz = C_nnz;
    C->sort();

    if (duplicates)
    {
        C->remove_duplicates();
    }
}


void CSRMatrix::add(CSRMatrix* B, CSRMatrix* C, bool duplicates)
{
    add(B, C, 1.0, duplicates);
}
void CSRMatrix::subtract(CSRMatrix* B, CSRMatrix* C)
{
    add(B, C, -1.0, true);
}




void add(BSRMatrix* B, BSRMatrix* C, double alpha, bool duplicates)
{
    int start, end;
    int C_nnz = nnz + B->nnz;
    C->cols.resize(C_nnz);
    C->data.resize(C_nnz);

    C_nnz = 0;
    C->rowptr[0] = 0;
    for (int i = 0; i < n_rows; i++)
    {
        start = rowptr[i];
        end = rowptr[i+1];
        std::copy(cols.begin() + start,
                cols.begin() + end,
                C->cols.begin() + C_nnz);
        if (data.size())
        {
            std::copy(data.begin() + start,
                    data.begin() + end,
                    C->data.begin() + C_nnz);
        }
        C_nnz += (end - start);

        start = B->rowptr[i];
        end = B->rowptr[i+1];
        std::copy(B->cols.begin() + start,
                B->cols.begin() + end,
                C->cols.begin() + C_nnz);
        if (B->data.size())
        {
            for (int j = 0; j < (end - start); j++)
            {
                std::array<double, b_size> block;
                for (int k = 0; k < b_size; k++)
                {
                    block[k] = alpha*B->data[start+j][k];
                }
                C->data[C_nnz+j] = block;
            }
        }
        C_nnz += (end - start);

        C->rowptr[i+1] = C_nnz;
    }
    C->nnz = C_nnz;
    C->sort();

    if (duplicates)
    {
        C->remove_duplicates();
    }
}


void BSRMatrix::add(BSRMatrix* B, BSRMatrix* C, bool duplicates)
{
    add(B, C, 1.0, duplicates);
}
void BSRMatrix::subtract(BSRMatrix* B, BSRMatrix* C)
{
    add(B, C, -1.0, true);
}


