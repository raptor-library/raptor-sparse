#include "src/core/matrix.hpp"

namespace raptor
{
    COOMatrix::COOMatrix(COOMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        rows.resize(A->nnz);
        std::copy(A->rows.begin(), A->rows.end(), rows.begin());

        cols.resize(A->nnz);
        std::copy(A->cols.begin(), A->cols.end(), cols.begin());

        if (A->data.size())
        {
            data.resize(A->nnz);
            std::copy(A->data.begin(), A->data.end(), data.begin());
        }
    }

    COOMatrix::COOMatrix(CSRMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    { 
        for (int i = 0; i < A->n_rows; i++)
        {
            int row_start = A->rowptr[i];
            int row_end = A->rowptr[i+1];
            for (int j = row_start; j < row_end; j++)
            {
                rows.push_back(i);
                cols.push_back(A->cols[j]);
                if (A->data.size())
                {
                    data.push_back(A->data[j]);
                }
            }
        }
    }

    COOMatrix::COOMatrix(CSCMatrix* A): Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        for (int i = 0; i < A->n_cols; i++)
        {
            int col_start = A->colptr[i];
            int col_end = A->colptr[i+1];
            for (int j = col_start; j < col_end; j++)
            {
                rows.push_back(A->rows[j]);
                cols.push_back(i);
                if (A->data.size())
                {
                    data.push_back(A->data[j]);
                }
            }
        }
    }

    COOMatrix* COOMatrix::transpose()
    {
        COOMatrix* T = new COOMatrix(n_cols, n_rows, cols, rows, data);
        return T;
    }

}
