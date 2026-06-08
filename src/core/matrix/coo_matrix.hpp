// Copyright (c) 2015-2026, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause
#ifndef RAPTOR_CORE_COO_MATRIX_HPP
#define RAPTOR_CORE_COO_MATRIX_HPP

#include "src/core/matrix.hpp"

/**************************************************************
 *****   COOMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in COO format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (COO)
 ***** sort()
 *****    Sorts the matrix by row, and by column within each row.
 ***** add_value(int row, int col, double val)
 *****     Adds val to position (row, col)
 ***** rows()
 *****     Returns std::vector<int>& containing the rows corresponding
 *****     to each nonzero
 ***** cols()
 *****     Returns std::vector<int>& containing the cols corresponding
 *****     to each nonzero
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 **************************************************************/
namespace raptor
{
  class CSRMatrix;
  class CSCMatrix;

  class COOMatrix : public Matrix
  {
    public:
        COOMatrix(int _nrows, int _ncols, int nnz_per_row = 1) :
                Matrix(_nrows, _ncols)
        {
            int _nnz = nnz_per_row * _nrows;
            if (_nnz)
            {
                rows.reserve(_nnz);
                cols.reserve(_nnz);
                data.reserve(_nnz);
            }
        }

        COOMatrix(int _nrows, int _ncols, double* _data) :
                Matrix(_nrows, _ncols)
        {
            for (int i = 0; i < n_rows; i++)
            {
                for (int j = 0; j < n_cols; j++)
                {
                    int pos = i*n_cols+j;
                    if (fabs(_data[pos]) > zero_tol)
                    {
                        rows.push_back(i);
                        cols.push_back(j);
                        data.push_back(_data[pos]);
                    }
                }
            }
        }

        COOMatrix(int _nrows, int _ncols, std::vector<int>& _rows,
                std::vector<int>& _cols, std::vector<double>& _data) :
            Matrix(_nrows, _ncols, _rows.size())
        {
            rows.resize(nnz);
            std::copy(_rows.begin(), _rows.end(), rows.begin());

            cols.resize(nnz);
            std::copy(_cols.begin(), _cols.end(), cols.begin());

            if (_data.size())
            {
                data.resize(nnz);
                std::copy(_data.begin(), _data.end(), data.begin());
            }
        }

        COOMatrix(COOMatrix* A);
        COOMatrix(CSRMatrix* A);
        COOMatrix(CSCMatrix* A);

        ~COOMatrix(){}

        COOMatrix* transpose();

        format_t format()
        {
            return COO;
        }

        template <typename T>
        void sort_helper(std::vector<T>& vals)
        {
            if (sorted || nnz == 0)
            {
                return;
            }

            vec_sort(rows, cols, vals);
            sorted = true;
            diag_first = false;
        }
        void sort()
        {
            sort_helper(data);
        }


        template <typename T>
        void move_diag_helper(std::vector<T>& vals)
        {
            if (diag_first || nnz == 0)
            {
                return;
            }

            int row_start, prev_row;
            int row, col;


            if (!sorted)
            {
                sort();
            }

            row_start = 0;
            prev_row = 0;
            for (int i = 0; i < nnz; i++)
            {
                row = rows[i];
                col = cols[i];
                if (row != prev_row)
                {
                    prev_row = row;
                    row_start = i;
                }
                else if (row == col)
                {
                    auto tmp = vals[i];
                    for (int j = i; j > row_start; j--)
                    {
                        cols[j] = cols[j-1];
                        vals[j] = vals[j-1];
                    }
                    cols[row_start] = row;
                    vals[row_start] = tmp;
                }
            }

            diag_first = true;
        }
        void move_diag()
        {
            move_diag_helper(data);
        }

        template <typename T>
        void remove_duplicates_helper(std::vector<T>& vals)
        {
            if (nnz == 0)
            {
                return;
            }

            if (!sorted)
            {
                sort();
            }

            int prev_row, prev_col, ctr;
            int row, col;

            // Remove duplicates (sum together)
            prev_row = rows[0];
            prev_col = cols[0];
            ctr = 1;
            for (int i = 1; i < nnz; i++)
            {
                row = rows[i];
                col = cols[i];
                if (row == prev_row && col == prev_col)
                {
                    append_vals(&vals[ctr - 1], &vals[i], b_size);
                }
                else
                { 
                    if (ctr != i)
                    {
                        rows[ctr] = row;
                        cols[ctr] = col;
                        vals[ctr] = vals[i];
                    }
                    ctr++;

                    prev_row = row;
                    prev_col = col;
                }
            }

            nnz = ctr;
        }
        void remove_duplicates()
        {
            remove_duplicates_helper(data);
        }

        void print()
        {
            for (int i = 0; i < nnz; i++)
            {
                printf("A[%d][%d] = %e\n", rows[i], cols[i], data[i]);
            }
        }

        void add_value(int row, int col, double val)
        {
            rows.push_back(row);
            cols.push_back(col);
            data.push_back(val);
            nnz++;
        }

        // Linear Algebra Methods
        void spmv(const double alpha, const double* x, 
                const double beta, double* b);
        void spmv_T(const double alpha, const double* x, 
                const double beta, double* b);

        std::vector<int> rows;
        std::vector<int> cols;
        std::vector<double> data;
  };
}

#endif
