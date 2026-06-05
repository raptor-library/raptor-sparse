#ifndef RAPTOR_CORE_BSC_MATRIX_HPP
#define RAPTOR_CORE_BSC_MATRIX_HPP

#include "src/core/matrix.hpp"

// TODO: Currently using column-major layout (e.g. fortran layout) 
// within blocks.  Is this correct, or do we want blocks to still 
// be row-major layout?

/**************************************************************
 *****   BSRMatrix Class (Inherits from CSRMatrix Base Class)
 **************************************************************
 ***** This class constructs a sparse matrix in BSR format.
 *****
 ***** Methods
 ***** -------
 ***** format() 
 *****    Returns the format of the sparse matrix (BSR)
 ***** sort()
 *****    Sorts the matrix.  Already in row-wise order, but sorts
 *****    the columns in each row.
 **************************************************************/
namespace raptor
{
  class CSCMatrix;
  class BSRMatrix;


  class BSCMatrix : public CSCMatrix
  {
    public:
        BSCMatrix(int _nrows, int _ncols, int _brows, 
                int _bcols, int _nnz = 0) 
        {
            n_rows = _nrows;
            n_cols = _ncols;
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;

            nnz = 0;
            sorted = false;
            diag_first = false;

            colptr.resize(n_cols + 1);
            if (_nnz)
            {
                rows.reserve(_nnz);
                data.reserve(_nnz);
            }
        }

        BSCMatrix(int _nrows, int _ncols, int _brows,
                int _bcols, double* _data)
        {
            n_rows = _nrows;
            n_cols = _ncols;
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;

            sorted = false;
            diag_first = false;

            colptr.resize(n_cols+1);

            nnz = 0;
            colptr[0] = 0;
            for (int i = 0; i < n_cols; i++)
            {
                for (int j = 0; j < n_rows; j++)
                {
                    double sum = 0;
                    std::array<double, b_size> block;
                    for (int bi = 0; bi < b_cols; bi++)
                    {
                        int col = i*b_cols+bi;
                        for (int bj = 0; bj < b_rows; bj++)
                        {
                            int row = j*b_rows+bj;
                            double val = _data[row*n_cols*b_cols+col];
                            block[bi*b_rows+bj] = val;
                            sum += fabs(val);
                        }
                    }
                    if (sum > zero_tol)
                    {
                        rows.push_back(j);
                        data.push_back(block);
                        nnz++;
                    }
                }
                colptr[i+1] = nnz;
            }
        }

        BSCMatrix(int _nrows, int _ncols, int _brows, int _bcols,
                std::vector<int>& _colptr, std::vector<int>& _rows,
                std::vector<double*> _data)
        {
            n_rows = _nrows;
            n_cols = _ncols;
            nnz = _rows.size();
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;

            colptr.resize(n_cols+1);
            std::copy(_colptr.begin(), _colptr.end(), colptr.begin());

            rows.resize(nnz);
            std::copy(_rows.begin(), _rows.end(), rows.begin());

            if (_data.size())
            {
                data.resize(nnz);
                for (int i = 0; i < nnz; i++)
                {
                    double* block = _data[i];
                    double* new_block = new double[b_size];
                    memcpy(new_block, block, b_size*sizeof(double));
                    data[i] = new_block;
                }
            }
        }

        BSCMatrix(CSCMatrix* A, int _brows, _bcols)
        {
            b_rows = _brows;
            b_cols = _bcols;
            b_size = b_rows * b_cols;
            n_rows = A->n_rows / b_rows;
            n_cols = A->n_cols / b_cols;
            colptr.resize(n_cols+1);

            // Map for each block col to accumulate block rows
            std::map<int, double*> block_map;

            colptr[0] = 0;
            for (int i = 0; i < A->n_cols; i++)
            {
                int bi = i / b_cols; 
                int bcol = i % b_cols;
                int col_start = A->colptr[i];
                int col_end = A->colptr[i+1];

                for (int j = col_start; j < col_end; j++)
                {
                    int row = A->rows[j];
                    int bj = row / b_rows;
                    int brow = row % b_rows;

                    // map key `bj` is block column to be stored in BSC matrix
                    if (block_map.find(bj) == block_map.end())
                        block_map[bj] = new double[b_size]();

                    // Add value into appropriate block / positon
                    block_map[bj][bcol*b_rows + brow] = A->data[j];
                }

                // Last col in a block
                if ((i+1)%b_cols == 0)
                {
                    for (auto& [bj, block] : block_map)
                    {
                        // Push cols and data
                        rows.push_back(bj);
                        data.push_back(block);
                    }
                    colptr[bi+1] = rows.size();
                    block_map.clear();
                }
            }
            nnz = rows.size();
        }


        ~BSCMatrix()
        {
            for (int i = 0; i < data.size(); i++)
            {
                delete data[i];
            }
        }


        BSCMatrix<b_cols, b_rows>* transpose()
        {
            BSRMatrix* T_bsr = new BSRMatrix(n_cols, n_rows, b_cols,
                    b_rows, colptr, rows, data);
            BSCMatrix* T = new BSCMatrix(T_bsr);
            delete T_bsr;
            return T;
        }


        format_t format()
        {
            return BSC;
        }

        void sort()
        {
            sort_helper(data);
        }

        void move_diag()
        {
            move_diag_helper(data);
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
                    double* = data[j];
                    for (int bcol = 0; bcol < b_cols; bcol++)
                    {
                        for (int brow = 0; brow < b_rows; brow++)
                        {
                            printf("A[%d][%d], Block [%d][%d] = %e\n", rows[j], col, brow, bcol, 
                                    val[bcol*b_rows+brow]);
                        }
                    }
                }
            }
        }

        // Underlying storage will be contiguous
        static constexpr int b_size = b_rows*b_cols;
        std::vector<std::array<double, b_size>> data;
  };
}

#endif


