#ifndef RAPTOR_CORE_CSC_MATRIX_HPP
#define RAPTOR_CORE_CSC_MATRIX_HPP

#include "src/core/matrix.hpp"

/**************************************************************
 *****   CSCMatrix Class (Inherits from Matrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in CSC format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (CSC)
 ***** sort()
 *****    Sorts the matrix.  Already in col-wise order, but sorts
 *****    the rows in each column.
 ***** add_value(int row, int col, double val)
 *****     TODO -- add this functionality
 ***** indptr()
 *****     Returns std::vector<int>& column pointer.  The ith element points to
 *****     the index of indices() corresponding to the first row to lie on 
 *****     column i.
 ***** indices()
 *****     Returns std::vector<int>& containing the rows corresponding
 *****     to each nonzero
 ***** data()
 *****     Returns std::vector<double>& containing the nonzero values
 **************************************************************/
namespace raptor
{
  class COOMatrix;
  class CSRMatrix;
  class BSCMatrix;

    class CSCMatrix : public Matrix
    {
        public:

        CSCMatrix(int _nrows, int _ncols, int _nnz = 0) : Matrix(_nrows, _ncols)
        {
            colptr.resize(n_cols+1);
            if (_nnz)
            {
                rows.reserve(_nnz);
                data.reserve(_nnz);
            }
        }

        CSCMatrix(int _nrows, int _ncols, double* _data) : Matrix(_nrows, _ncols)
        {
            colptr.resize(n_cols+1);

            nnz = 0;
            colptr[0] = 0;
            for (int i = 0; i < n_cols; i++)
            {
                for (int j = 0; j < n_rows; j++)
                {
                    int pos = i*n_cols+j;
                    if (fabs(_data[pos]) > zero_tol)
                    {
                        rows.push_back(j);
                        data.push_back(_data[pos]);
                        nnz++;
                    }
                }
                colptr[i+1] = nnz;
            }
        }

        CSCMatrix(int _nrows, int _ncols, std::vector<int>& _colptr,
                std::vector<int>& _rows, std::vector<double>& _data) :
            Matrix(_nrows, _ncols, _rows.size())
        {
            colptr.resize(n_cols+1);
            std::copy(_colptr.begin(), _colptr.end(), colptr.begin());

            rows.resize(nnz);
            std::copy(_rows.begin(), _rows.end(), rows.begin());

            if (_data.size())
            {
                data.resize(nnz);
                std::copy(_data.begin(), _data.end(), data.begin());
            }
        }

        CSCMatrix(CSCMatrix* A);
        CSCMatrix(COOMatrix* A);
        CSCMatrix(CSRMatrix* A);
        
        ~CSCMatrix() {}

        CSCMatrix* transpose();

        format_t format()
        {
            return CSC;
        }

        template <typename T>
        void sort_helper(std::vector<T>& vals)
        {
            if (sorted || nnz == 0)
            {
                return;
            }

            int start, end, size;
    
            // Sort the columns of each col (and data accordingly) and remove
            // duplicates (summing values together)
            for (int col = 0; col < n_cols; col++)
            {
                start = colptr[col];
                end = colptr[col+1];
                size = end - start;
                if (size == 0)
                {
                    continue;
                }

                if (vals.size())
                    vec_sort(rows, vals, start, end);
                else
                    std::sort(rows.begin() + start, rows.begin() + end);
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

            int start, end;

            if (!sorted)
            {
                sort();
            }

            T tmp;
            for (int i = 0; i < n_cols; i++)
            {
                start = colptr[i];
                end = colptr[i+1];
                for (int j = start; j < end; j++)
                {
                    if (vals.size())
                    {
                        tmp = vals[j];
                    }
                    for (int k = j; k > start; k--)
                    {
                        rows[k] = rows[k-1];
                        if (vals.size())
                        {
                            vals[k] = vals[k-1];
                        }
                    }
                    rows[start] = i;
                    if (vals.size())
                    {
                        vals[start] = tmp;
                    }
                    break;
                }
            }
        }
        void move_diag()
        {
            move_diag_helper(data);
        }

        template <typename T>
        void remove_duplicates_helper(std::vector<T>& vals)
        {
            int orig_start, orig_end;
            int new_start;
            int row, prev_row;
            int ctr, col_size;

            if (!sorted)
            {
                sort();
            }

            orig_start = colptr[0];
            for (int col = 0; col < n_cols; col++)
            {
                new_start = colptr[col];
                orig_end = colptr[col+1];
                col_size = orig_end - orig_start;
                if (col_size == 0) 
                {
                    orig_start = orig_end;
                    colptr[col+1] = colptr[col];
                    continue;
                }

                // Remove Duplicates
                row = rows[orig_start];
                rows[new_start] = row;
                vals[new_start] = vals[orig_start];
                prev_row = row;
                ctr = 1;
                for (int j = orig_start + 1; j < orig_end; j++)
                {
                    row = rows[j];
                    if (row == prev_row)
                    {
                        append_vals(&vals[ctr - 1 + new_start], &vals[j], b_size);
                    }
                    else
                    {
                        if (abs_val(vals[ctr - 1 + new_start], b_size) < zero_tol)
                        {
                            ctr--;
                        }

                        rows[ctr + new_start] = row;
                        vals[ctr + new_start] = vals[j];
                        ctr++;
                        prev_row = row;
                    }
                }
                if (abs_val(vals[ctr - 1 + new_start], b_size) < zero_tol)
                {
                    ctr--;
                }

                orig_start = orig_end;
                colptr[col+1] = colptr[col] + ctr;
            }
            nnz = colptr[n_cols];
            rows.resize(nnz);
            vals.resize(nnz);
        }
        void remove_duplicates()
        {
            remove_duplicates_helper(data);
        }

        void print()
        {
            int col_start, col_end;
            for (int col = 0; col < n_cols; col++)
            {
                col_start = colptr[col];
                col_end = colptr[col+1];
                for (int j = col_start; j < col_end; j++)
                {
                    printf("A[%d][%d] = %e\n", rows[j], col, data[j]);                
                }
            }
        }

        // Linear Algebra Methods
        void spmv(const double alpha, const double* x, 
                const double beta, double* b);

        std::vector<int> colptr;
        std::vector<int> rows;
        std::vector<double> data;
    };
}

#endif
