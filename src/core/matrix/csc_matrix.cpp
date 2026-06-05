#include "src/core/matrix.hpp"

namespace raptor
{
    CSCMatrix::CSCMatrix(CSCMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        colptr.resize(n_cols+1);
        std::copy(A->colptr.begin(), A->colptr.end(), colptr.begin());

        rows.resize(nnz);
        std::copy(A->rows.begin(), A->rows.end(), rows.begin());

        if (A->data.size())
        {
            data.resize(nnz);
            std::copy(A->data.begin(), A->data.end(), data.begin());
        }

    }

    CSCMatrix::CSCMatrix(COOMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        colptr.resize(n_cols+1);
        rows.resize(nnz);

        if (A->data.size())
        {
            data.resize(nnz);
        }

        std::fill(colptr.begin(), colptr.end(), 0);
        for (int i = 0; i < nnz; i++)
        {
            colptr[A->cols[i]+1]++;
        }
        for (int i = 0; i < n_cols; i++)
        {
            colptr[i+1] += colptr[i];
        }

        std::vector<int> ctr(n_cols, 0);
        for (int i  = 0; i < nnz; i++)
        {
            int row = A->rows[i];
            int col = A->cols[i];
            int index = colptr[col] + ctr[col]++;
            rows[index] = row;
            if (A->data.size())
            {
                data[index] = A->data[i];
            }
        }
    }

    CSCMatrix::CSCMatrix(CSRMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        colptr.resize(n_cols+1);
        rows.resize(nnz);
        if (A->data.size())
        {
            data.resize(nnz);
        }

        std::fill(colptr.begin(), colptr.end(), 0);
        for (int i = 0; i < nnz; i++)
        {
            colptr[A->cols[i]+1]++;
        }
        for (int i = 0; i < n_cols; i++)
        {
            colptr[i+1] += colptr[i];
        }

        std::vector<int> ctr(n_cols, 0);
        for (int i  = 0; i < n_rows; i++)
        {
            int row_start = A->rowptr[i];
            int row_end = A->rowptr[i+1];
            for (int j = row_start; j < row_end; j++)
            {
                int col = A->cols[j];
                int index = colptr[col] + ctr[col]++;
                rows[index] = i;
                if (A->data.size())
                {
                    data[index] = A->data[j];
                }
            }
        }
    }

    CSCMatrix* CSCMatrix::transpose()
    {
        CSRMatrix* T_csr = new CSRMatrix(n_cols, n_rows, colptr, rows, data);
        CSCMatrix* T = new CSCMatrix(T_csr);
        delete T_csr;
        return T;
    }
}
