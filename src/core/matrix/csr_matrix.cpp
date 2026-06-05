#include "src/core/matrix.hpp"

namespace raptor
{

    CSRMatrix::CSRMatrix(CSRMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        rowptr.resize(n_rows+1);
        std::copy(A->rowptr.begin(), A->rowptr.end(), rowptr.begin());

        cols.resize(nnz);
        std::copy(A->cols.begin(), A->cols.end(), cols.begin());

        if (A->data.size())
        {
            data.resize(nnz);
            std::copy(A->data.begin(), A->data.end(), data.begin());
        }
    }

    CSRMatrix::CSRMatrix(COOMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        rowptr.resize(n_rows+1);
        cols.resize(nnz);

        if (A->data.size())
        {
            data.resize(nnz);
        }

        std::fill(rowptr.begin(), rowptr.end(), 0);
        for (int i = 0; i < nnz; i++)
        {
            rowptr[A->rows[i]+1]++;
        }
        for (int i = 0; i < n_rows; i++)
        {
            rowptr[i+1] += rowptr[i];
        }

        std::vector<int> ctr(n_rows, 0);
        for (int i = 0; i < nnz; i++)
        {
            int row = A->rows[i];
            int col = A->cols[i];
            int index = rowptr[row] + ctr[row]++;
            cols[index] = col;
            if (A->data.size())
            {
                data[index] = A->data[i];
            }
        }
    }

    CSRMatrix::CSRMatrix(CSCMatrix* A) : Matrix(A->n_rows, A->n_cols, A->nnz)
    {
        rowptr.resize(n_rows+1);
        cols.resize(nnz);
        if (A->data.size())
        {
            data.resize(nnz);
        }

        std::fill(rowptr.begin(), rowptr.end(), 0);
        for (int i = 0; i < nnz; i++)
        {
            rowptr[A->rows[i]+1]++;
        }
        for (int i = 0; i < n_rows; i++)
        {
            rowptr[i+1] += rowptr[i];
        }

        std::vector<int> ctr(n_rows, 0);
        for (int i = 0; i < n_rows; i++)
        {
            int col_start = A->colptr[i];
            int col_end = A->colptr[i+1];
            for (int j = col_start; j < col_end; j++)
            {
                int row = A->rows[j];
                int index = rowptr[row] + ctr[row]++;
                cols[index] = i;
                if (A->data.size())
                {
                    data[index] = A->data[j];
                }
            }
        }
    }

    CSRMatrix::CSRMatrix(BSRMatrix* A)
    {
        b_rows = 1;
        b_cols = 1;
        b_size = 1;

        n_rows = A->n_rows * A->b_rows;
        n_cols = A->n_cols * A->b_cols;

        rowptr.resize(n_rows+1);
        rowptr[0] = 0;
        for (int i = 0; i < A->n_rows; i++)
        {
            int row_start = A->rowptr[i];
            int row_end = A->rowptr[i+1];
            for (int br = 0; br < A->b_rows; br++)
            {
                for (int j = row_start; j < row_end; j++)
                {
                    for (int bc = 0; bc < A->b_cols; bc++)
                    {
                        double val = A->data[j][br*A->b_cols + bc];
                        if (fabs(val) > zero_tol)
                        {
                            int col = A->cols[j];
                            cols.push_back(col*A->b_cols + bc);
                            data.push_back(val);
                        }
                    }
                }
                rowptr[i*A->b_rows + br+1] = cols.size();
            }
        }
        nnz = cols.size();
        sorted = false;
        diag_first = false;
    }

    CSRMatrix* CSRMatrix::transpose()
    {
        CSCMatrix* T_csc = new CSCMatrix(n_cols, n_rows, rowptr, cols, data);
        CSRMatrix* T = new CSRMatrix(T_csc);
        delete T_csc;
        return T;
    }


}
