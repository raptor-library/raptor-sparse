#ifndef RAPTOR_CORE_CSR_MATRIX_HPP
#define RAPTOR_CORE_CSR_MATRIX_HPP

#include "src/core/matrix/matrix_base.hpp"

/**************************************************************
 *****   CSRMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in CSR format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (CSR)
 ***** sort()
 *****    Sorts the matrix.  Already in row-wise order, but sorts
 *****    the columns in each row.
 ***** add_value(int row, int col, double val)
 *****     TODO -- add this functionality
 ***** indptr()
 *****     Returns std::vector<int>& row pointer.  The ith element points to
 *****     the index of indices() corresponding to the first column to lie on 
 *****     row i.
 ***** indices()
 *****     Returns std::vector<int>& containing the cols corresponding
 *****     to each nonzero
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 **************************************************************/
namespace raptor
{
  class COOMatrix;
  class CSCMatrix;
  class BSRMatrix;

  class CSRMatrix : public Matrix
  {
    public:
        CSRMatrix(int _nrows, int _ncols, int _nnz = 0): Matrix(_nrows, _ncols)
        {
            rowptr.resize(n_rows + 1);
            if (_nnz)
            {
                cols.reserve(_nnz);
                data.reserve(_nnz);
            }
        }

        CSRMatrix(int _nrows, int _ncols, double* _data) :
                Matrix(_nrows, _ncols)
        {
            rowptr.resize(n_rows+1);

            nnz = 0;
            rowptr[0] = 0;
            for (int i = 0; i < n_rows; i++)
            {
                for (int j = 0; j < n_cols; j++)
                {
                    int pos = i*n_cols+j;
                    if (fabs(_data[pos]) > zero_tol)
                    {
                        cols.push_back(j);
                        data.push_back(_data[pos]);
                        nnz++;
                    }
                }
                rowptr[i+1] = nnz;
            }
        }

        CSRMatrix(int _nrows, int _ncols, std::vector<int>& _rowptr,
                std::vector<int>& _cols, std::vector<double>& _data) :
            Matrix(_nrows, _ncols, _cols.size())
        {
            rowptr.resize(n_rows+1);
            std::copy(_rowptr.begin(), _rowptr.end(), rowptr.begin());

            cols.resize(nnz);
            std::copy(_cols.begin(), _cols.end(), cols.begin());

            if (_data.size())
            {
                data.resize(nnz);
                std::copy(_data.begin(), _data.end(), data.begin());
            }
        }

        CSRMatrix(CSRMatrix* A);
        CSRMatrix(COOMatrix* A);
        CSRMatrix(CSCMatrix* A);

        ~CSRMatrix(){}

        CSRMatrix* transpose();

        format_t format()
        {
            return CSR;
        }

        template <typename T>
        void sort_helper(std::vector<T>& vals)
        {
            if (sorted || nnz == 0)
            {
                return;
            }

            int start, end, size;
            for (int row = 0; row < n_rows; row++)
            {
                start = rowptr[row];
                end = rowptr[row+1];
                size = end - start;
                if (size == 0) 
                {
                    continue;
                }
                if (vals.size())
                    vec_sort(cols, vals, start, end);
                else
                    std::sort(cols.begin() + start, cols.begin() + end);
            }

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

            int start, end, col;
            T tmp;
            // Move diagonal values to beginning of each row
            for (int i = 0; i < n_rows; i++)
            {
                start = rowptr[i];
                end = rowptr[i+1];
                for (int j = start; j < end; j++)
                {
                    col = cols[j];
                    if (col == i)
                    {
                        if (vals.size())
                        { 
                            tmp = vals[j];
                        }
                        for (int k = j; k > start; k--)
                        {
                            cols[k] = cols[k-1];
                            if (vals.size())
                            {
                                vals[k] = vals[k-1];
                            }
                        }
                        cols[start] = i;
                        if (vals.size())
                        {
                            vals[start] = tmp;
                        }
                        break;
                    }
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
            if (!sorted)
            {
                sort();
            }

            int orig_start, orig_end, new_start;
            int col, prev_col, ctr, row_size;

            orig_start = rowptr[0];
            for (int row = 0; row < n_rows; row++)
            {
                new_start = rowptr[row];
                orig_end = rowptr[row+1];
                row_size = orig_end - orig_start;
                if (row_size == 0) 
                {
                    orig_start = orig_end;
                    rowptr[row+1] = rowptr[row];
                    continue;
                }

                // Remove Duplicates
                col = cols[orig_start];
                cols[new_start] = col;
                if (vals.size())
                {
                    vals[new_start] = vals[orig_start];
                }
                prev_col = col;
                ctr = 1;
                for (int j = orig_start + 1; j < orig_end; j++)
                {
                    col = cols[j];
                    if (col == prev_col)
                    {
                        if (vals.size())
                        {
                            append_vals(&vals[ctr - 1 + new_start], &vals[j], b_size);
                        }
                    }
                    else
                    {
                        if (vals.size())
                        {
                            if (abs_val(vals[ctr - 1 + new_start], b_size) < zero_tol)
                            {
                                ctr--;
                            }
                            vals[ctr + new_start] = vals[j];
                        }

                        cols[ctr + new_start] = col;
                        ctr++;
                        prev_col = col;
                    }
                }
                if (vals.size() && abs_val(vals[ctr - 1 + new_start], b_size) < zero_tol)
                {
                    ctr--;
                }

                orig_start = orig_end;
                rowptr[row+1] = rowptr[row] + ctr;
            }
            nnz = rowptr[n_rows];
            cols.resize(nnz);
            if (vals.size())
            {
                vals.resize(nnz);
            }
        }
        void remove_duplicates()
        {
            remove_duplicates_helper(data);
        }

        void print()
        {
            int row_start, row_end;
            for (int row = 0; row < n_rows; row++)
            {
                row_start = rowptr[row];
                row_end = rowptr[row+1];
                for (int j = row_start; j < row_end; j++)
                {
                    printf("A[%d][%d] = %e\n", row, cols[j], data[j]);
                }
            }
        }


        // Linear Algebra Methods
        void spmv(const double alpha, const double* x, 
                const double beta, double* b);
        void add(CSRMatrix* B, CSRMatrix* C, bool duplicates);
        void subtract(CSRMatrix* B, CSRMatrix* C);

        std::vector<int> rowptr;
        std::vector<int> cols;
        std::vector<double> data;
  };
}

#endif
