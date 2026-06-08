// Copyright (c) 2015-2017, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause

#include "src/core/matrix.hpp"

using namespace raptor;

void CSRMatrix::spmv(const double alpha, const double* x, 
        const double beta, double* b)
{
    int start, end;
    double val;
    for (int i = 0; i < n_rows; i++)
    {
        start = rowptr[i];
        end = rowptr[i+1];
        val = 0;
        for (int j = start; j < end; j++)
        {
            val += alpha * data[j] * x[cols[j]];
        }
        b[i] = val + beta*b[i];
    }
}
void CSRMatrix::spmv_T(const double alpha, const double* x, 
        const double beta, double* b)
{
    int start, end;
    for (int i = 0; i < n_cols; i++)
    {
        b[i] = beta*b[i];
    }

    for (int i = 0; i < n_rows; i++)
    {
        start = rowptr[i];
        end = rowptr[i+1];
        for (int j = start; j < end; j++)
        {
            b[cols[j]] += alpha * data[j] * x[i];
        }
    }
}
void CSCMatrix::spmv(const double alpha, const double* x, 
        const double beta, double* b)
{
    int start, end;

    for (int i = 0; i < n_rows; i++)
    {
        b[i] = beta*b[i];
    }
    for (int i = 0; i < n_cols; i++)
    {
        start = colptr[i];
        end = colptr[i+1];
        for (int j = start; j < end; j++)
        {
            b[rows[j]] += alpha * data[j] * x[i];
        }
    }
}
void CSCMatrix::spmv_T(const double alpha, const double* x, 
        const double beta, double* b)
{
    int start, end;
    double val;

    for (int i = 0; i < n_cols; i++)
    {
        start = colptr[i];
        end = colptr[i+1];
        val = 0;
        for (int j = start; j < end; j++)
        {
            val += alpha * data[j] * x[rows[j]];
        }
        b[i] = val + beta*b[i];
    }
}
void COOMatrix::spmv(const double alpha, const double* x,
        const double beta, double* b)
{
    for (int i = 0; i < n_rows; i++)
    {
        b[i] = beta*b[i];
    }
    for (int i = 0; i < nnz; i++)
    {
        b[rows[i]] += alpha * data[i] * x[cols[i]];
    }
}
void COOMatrix::spmv_T(const double alpha, const double* x,
        const double beta, double* b)
{
    for (int i = 0; i < n_cols; i++)
    {
        b[i] = beta*b[i];
    }
    for (int i = 0; i < nnz; i++)
    {
        b[cols[i]] += alpha * data[i] * x[rows[i]];
    }
}
void BSRMatrix::spmv(const double alpha, const double* x,
        const double beta, double* b)
{
    int start, end, first_row, first_col, idx;
    double val;
    for (int i = 0; i < n_rows; i++)
    {
        start = rowptr[i];
        end = rowptr[i+1];
        first_row = i*b_rows;
        for (int bi = 0; bi < b_rows; bi++)
        {
            val = 0;
            idx = bi * b_cols;
            for (int j = start; j < end; j++)
            {
                first_col = cols[j]*b_cols;
                double* block = data[j];
                for (int bj = 0; bj < b_cols; bj++)
                {
                    val += alpha * block[bi*b_cols+bj] * x[first_col+bj];
                }
            }
            b[first_row+bi] = val + beta*b[first_row+bi];
        }
    }
}
void BSRMatrix::spmv_T(const double alpha, const double* x,
        const double beta, double* b)
{
    int start, end, first_row, first_col, idx;
    for (int i = 0; i < n_cols*b_cols; i++)
    {
       b[i] = beta*b[i];
    } 
    for (int i = 0; i < n_rows; i++)
    {
        start = rowptr[i];
        end = rowptr[i+1];
        first_row = i*b_rows;
        for (int bi = 0; bi < b_rows; bi++)
        {
            idx = bi * b_cols;
            for (int j = start; j < end; j++)
            {
                first_col = cols[j]*b_cols;
                double* block = data[j];
                for (int bj = 0; bj < b_cols; bj++)
                {
                    b[first_col+bj] += alpha * block[bi*b_cols+bj] * x[first_row+bi];
                }
            }
        }
    }
}
void BSCMatrix::spmv(const double alpha, const double* x,
        const double beta, double* b)
{
    int start, end, first_row, first_col, idx;

    for (int i = 0; i < n_rows*b_rows; i++)
    {
        b[i] = beta * b[i];
    }
    for (int i = 0; i < n_cols; i++)
    {
        start = colptr[i];
        end = colptr[i+1];
        first_col = i*b_cols;
        for (int bi = 0; bi < b_cols; bi++)
        {
            idx = bi * b_rows;
            for (int j = start; j < end; j++)
            {
                first_row = rows[j]*b_rows;
                double* block = data[j];
                for (int bj = 0; bj < b_rows; bj++)
                {
                    b[first_row+bj] += alpha * block[bi*b_rows+bj] * x[first_col+bi];
                }
            }
        }
    }
}
void BSCMatrix::spmv_T(const double alpha, const double* x,
        const double beta, double* b)
{
    int start, end, first_row, first_col, idx;
    double val;

    for (int i = 0; i < n_cols; i++)
    {
        start = colptr[i];
        end = colptr[i+1];
        first_col = i*b_cols;
        for (int bi = 0; bi < b_cols; bi++)
        {
            val = 0;
            idx = bi * b_rows;
            for (int j = start; j < end; j++)
            {
                first_row = rows[j]*b_rows;
                double* block = data[j];
                for (int bj = 0; bj < b_rows; bj++)
                {
                    val+= alpha * block[bi*b_rows+bj] * x[first_row+bj];
                }
            }
            b[first_col+bi] = val + beta*b[first_col+bi];
        }
    }
}

void Matrix::mult(const double* x, double* b)
{
    spmv(1.0, x, 0.0, b);
}
void Matrix::mult_T(const double* x, double * b)
{
    spmv_T(1.0, x, 0.0, b);
}
void Matrix::residual(const double* x, const double* b, double* r)
{

    memcpy(r, b, n_rows * sizeof(double));
    spmv(1.0, x, -1.0, r);
}



