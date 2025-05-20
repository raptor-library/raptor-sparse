// Copyright (c) 2015-2017, RAPtor Developer Team
// License: Simplified BSD, http://opensource.org/licenses/BSD-2-Clause
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "raptor-sparse.hpp"

#ifdef USING_MFEM
#include "raptor/external/mfem_wrapper.hpp"
#endif

using namespace raptor;

#define BS 25


void mult_append(CSRMatrix* A, double* x, double* b)
{
    int start, end, col;
    double A_val;
    for (int i = 0; i < A->n_rows; i++)
    {
        start = A->idx1[i];
        end = A->idx1[i+1];
        for (int j = start; j < end; j++)
        {
            A_val = A->vals[j];
            col = A->idx2[j];
            for (int k = 0; k < BS; k++)
            {
                b[i*BS+k] += A_val * x[col*BS+k];
            }
        }
    }
}

void mult(CSRMatrix* A, double* x, double* b)
{
    for (int i = 0; i < A->n_rows; i++)
        for (int j = 0; j < BS; j++)
            b[i*BS+j] = 0;

    mult_append(A, x, b);
}



void mult_append(CSCMatrix* A, double* x, double* b)
{
    int start, end, row;
    double A_val;

    for (int col = 0; col < A->n_cols; col++)
    {
        start = A->idx1[col];
        end = A->idx1[col+1];
        for (int j = start; j < end; j++)
        {
            A_val = A->vals[j];
            row = A->idx2[j];
            for (int k = 0; k < BS; k++)
            {
                b[row*BS+k] += A_val * x[col*BS+k];
            }
        }
    }
}

void mult_append_msg(CSCMatrix* A, double* x, double* b, int start, int end)
{
    int col_start, col_end, row;
    double A_val;

	for (int col = start; col < end; col++)
	{
        col_start = A->idx1[col];
		col_end = A->idx1[col+1];
		for (int j = col_start; j < col_end; j++)
		{
			A_val = A->vals[j];
			row = A->idx2[j];
			for (int k = 0; k < BS; k++)
			{
				b[row*BS+k] += A_val * x[col*BS+k];
			}
		}
	}
}


void init_sends(ParCSRMatrix* A, double* values)
{
    int proc, start, end, idx;

    A->comm->send_data->buffer.resize(A->comm->send_data->size_msgs*BS);

    // Initialize Sends
    for (int i = 0; i < A->comm->send_data->num_msgs; i++)
    {
        proc = A->comm->send_data->procs[i];
        start = A->comm->send_data->indptr[i];
        end = A->comm->send_data->indptr[i+1];
        for (int j = start; j < end; j++)
        {
            idx = A->comm->send_data->indices[j];
            memcpy(&(A->comm->send_data->buffer[j*BS]), &(values[idx*BS]), BS);
        }
        MPI_Isend(&(A->comm->send_data->buffer[start*BS]), (end - start)*BS, MPI_DOUBLE, 
                proc, A->comm->key, A->comm->mpi_comm, &(A->comm->send_data->requests[i]));
   }
}

void init_recvs(ParCSRMatrix* A)
{
    int proc, start, end;

    A->comm->recv_data->buffer.resize(A->comm->recv_data->size_msgs*BS);

    for (int i = 0; i < A->comm->recv_data->num_msgs; i++)
    {
        proc = A->comm->recv_data->procs[i];
        start = A->comm->recv_data->indptr[i];
        end = A->comm->recv_data->indptr[i+1];
        MPI_Irecv(&(A->comm->recv_data->buffer[start*BS]), (end - start)*BS, MPI_DOUBLE, 
            proc, A->comm->key, A->comm->mpi_comm, &(A->comm->recv_data->requests[i]));
    }
}

void wait_sends(ParCSRMatrix* A)
{
    if (A->comm->send_data->num_msgs)
        MPI_Waitall(A->comm->send_data->num_msgs, 
                A->comm->send_data->requests.data(),
                MPI_STATUSES_IGNORE);
}

void wait_recvs(ParCSRMatrix* A)
{
    if (A->comm->recv_data->num_msgs)
        MPI_Waitall(A->comm->recv_data->num_msgs,
                A->comm->recv_data->requests.data(),
                MPI_STATUSES_IGNORE);
}

void standard_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;
    CSRMatrix* A_off = (CSRMatrix*)A->off_proc;

    init_recvs(A);
    init_sends(A, x.local.data());
    mult(A_on, x.local.data(), b.local.data());
    wait_recvs(A);
    mult_append(A_off, A->comm->recv_data->buffer.data(), b.local.data());
    wait_sends(A);
}

void early_comp_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;
    CSRMatrix* A_off = (CSRMatrix*)A->off_proc;

    init_recvs(A);
    init_sends(A, x.local.data());
    mult(A_on, x.local.data(), b.local.data());
    wait_recvs(A);
    mult_append(A_off, A->comm->recv_data->buffer.data(), b.local.data());
    wait_sends(A);
}

void CSC_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b, CSCMatrix* A_off_csc)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;

    init_recvs(A);
    init_sends(A, x.local.data());
    mult(A_on, x.local.data(), b.local.data());
    wait_recvs(A);
    mult_append(A_off_csc, A->comm->recv_data->buffer.data(), b.local.data());
    wait_sends(A);
}


void CSC_early_recv_spmv(ParCSRMatrix* A, ParVector& x, ParVector& b, CSCMatrix* A_off_csc)
{
    CSRMatrix* A_on = (CSRMatrix*)A->on_proc;

    int start, end, idx;

    MPI_Status status;
    int flag = 0;

    init_recvs(A);
    init_sends(A, x.local.data());
    mult(A_on, x.local.data(), b.local.data());

    int n_req = A->comm->recv_data->num_msgs;
    int next_n_req = 0;
    std::vector<int> req_idx(A->comm->recv_data->num_msgs);
    std::iota(req_idx.begin(), req_idx.end(), 0);
    while (n_req)
	{
        next_n_req = 0;
		for (int i = 0; i < n_req; i++)
        {
            idx = req_idx[i];
		    MPI_Test(&(A->comm->recv_data->requests[idx]), &flag, &status);
            if (flag)
			{
                start = A->comm->recv_data->indptr[idx];
				end = A->comm->recv_data->indptr[idx+1];
				mult_append_msg(A_off_csc, A->comm->recv_data->buffer.data(), b.local.data(), start, end);
			}
			else
			{
				req_idx[next_n_req++] = idx;
			}
		}
        n_req = next_n_req;
	}
    wait_sends(A);
}


using namespace raptor;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int dim=0;
    int n = 5;
    int system = 0;

    if (argc > 1)
    {
        system = atoi(argv[1]);
    }

    ParCSRMatrix* A=NULL;
    ParVector x;
    ParVector b;

    double t0, tfinal;

    int cache_len = 10000;

    std::vector<double> cache_array(cache_len);

    if (system < 2)
    {
        double* stencil = NULL;
        std::vector<int> grid;
        if (argc > 2)
        {
            n = atoi(argv[2]);
        }

        if (system == 0)
        {
            dim = 3;
            grid.resize(dim, n);
            stencil = laplace_stencil_27pt();
        }
        else if (system == 1)
        {
            dim = 2;
            grid.resize(dim, n);
            double eps = 0.001;
            double theta = M_PI/8.0;
            if (argc > 3)
            {
                eps = atof(argv[3]);
                if (argc > 4)
                {
                    theta = atof(argv[4]);
                }
            }
            stencil = diffusion_stencil_2d(eps, theta);
        }
        A = par_stencil_grid(stencil, grid.data(), dim);
        delete[] stencil;
    }
#ifdef USING_MFEM
    if (system < 2)
    {
        double* stencil = NULL;
        std::vector<int> grid;
        if (argc > 2)
        {
            n = atoi(argv[2]);
        }

        if (system == 0)
        {
            dim = 3;
            grid.resize(dim, n);
            stencil = laplace_stencil_27pt();
        }
        else if (system == 1)
        {
            dim = 2;
            grid.resize(dim, n);
            double eps = 0.001;
            double theta = M_PI/4.0;
            if (argc > 3)
            {
                eps = atof(argv[3]);
                if (argc > 4)
                {
                    theta = atof(argv[4]);
                }
            }
            stencil = diffusion_stencil_2d(eps, theta);
        }
        A = par_stencil_grid(stencil, grid.data(), dim);
        delete[] stencil;
    }
#endif
    else if (system == 3)
    {
        const char* file = "../../test_data/rss_A0.pm";
        if (argc > 2)
        {
            file = argv[2];
        }
        A = readParMatrix(file);
    }

    if (system != 2)
    {
        x = ParVector(A->global_num_cols*BS, A->on_proc_num_cols*BS);
        b = ParVector(A->global_num_rows*BS, A->local_num_rows*BS);
        x.set_const_value(1.0);
    }

    ParVector b_new = ParVector(A->global_num_rows*BS, A->local_num_rows*BS);
    CSCMatrix* A_off_csc = A->off_proc->to_CSC();

    /************************
    **** Correctness Checks
    ************************/
    b.set_const_value(0.0);
    standard_spmv(A, x, b);

    b_new.set_const_value(0.0);
    early_comp_spmv(A, x, b_new);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");
    
    b_new.set_const_value(0.0);
    CSC_spmv(A, x, b_new, A_off_csc);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");

    b_new.set_const_value(0.0);
    CSC_early_recv_spmv(A, x, b_new, A_off_csc);
    for (int i = 0; i < A->local_num_rows; i++)
        if (fabs(b[i] - b_new[i]) > 1e-06)
            printf("Difference in RESULTS!\n");

    int n_iters;
    MPI_Barrier(MPI_COMM_WORLD);

    /************************
    **** Time Standard SpMV
    ************************/
    // Warm-up
    standard_spmv(A, x, b);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    standard_spmv(A, x, b);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            standard_spmv(A, x, b);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("Standard SpMV Time: %e\n", t0);
    }

    /************************
    **** Time Early Computation SpMV
    ************************/
    // Warm-up
    early_comp_spmv(A, x, b_new);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    early_comp_spmv(A, x, b_new);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            early_comp_spmv(A, x, b_new);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("Early Comp SpMV Time: %e\n", t0);
    }

    /************************
    **** Time CSC SpMV
    ************************/
    // Warm-up
    CSC_spmv(A, x, b_new, A_off_csc);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    CSC_spmv(A, x, b_new, A_off_csc);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            CSC_spmv(A, x, b_new, A_off_csc);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("CSC SpMV Time: %e\n", t0);
    }


    /************************
    **** Time Early Recv CSC SpMV
    ************************/
    // Warm-up
    CSC_early_recv_spmv(A, x, b_new, A_off_csc);

    // Calculate iterations
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    CSC_early_recv_spmv(A, x, b_new, A_off_csc);
    tfinal = MPI_Wtime() - t0;
    MPI_Allreduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    n_iters = 1.0 / t0 + 1;

    for (int test = 0; test < 5; test++)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
        for (int i = 0; i < n_iters; i++)
        {
            CSC_early_recv_spmv(A, x, b_new, A_off_csc);
        }
        tfinal = (MPI_Wtime() - t0) / n_iters;
        MPI_Reduce(&tfinal, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("Early Recv CSC SpMV Time: %e\n", t0);
    }

    delete A_off_csc;
    delete A;

    MPI_Finalize();

    return 0;
}

