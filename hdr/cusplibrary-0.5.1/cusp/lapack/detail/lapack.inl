/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*! \file lapack.inl
 *  \brief Definition of lapack interface routines
 */

#include <cusp/array1d.h>
#include <cusp/array2d.h>

#include <cusp/lapack/detail/stubs.h>

namespace cusp
{
namespace lapack
{

template<typename Array2d, typename Array1d>
void getrf( Array2d& A, Array1d& piv )
{
    typedef typename Array2d::value_type ValueType;
    typedef typename Array1d::value_type IndexType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    if(piv.size() != A.num_cols)
        piv.resize(A.num_cols);

    lapack_int m    = A.num_cols;
    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    ValueType *a    = thrust::raw_pointer_cast(&A(0,0));
    IndexType *ipiv = thrust::raw_pointer_cast(&piv[0]);

    lapack_int info = detail::getrf(order, m, n, a, lda, ipiv);

    if( info != 0 )
        throw cusp::runtime_exception("getrf failed");
}

template<typename Array2d>
void potrf( Array2d& A, char uplo )
{
    typedef typename Array2d::value_type ValueType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    ValueType *a    = thrust::raw_pointer_cast(&A(0,0));
    lapack_int info = detail::potrf(order, uplo, n, a, lda);

    if( info != 0 )
        throw cusp::runtime_exception("potrf failed");
}

template<typename Array2d, typename Array1d>
void sytrf( Array2d& A, Array1d& piv, char uplo )
{
    typedef typename Array2d::value_type ValueType;
    typedef typename Array1d::value_type IndexType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    if(piv.size() != A.num_cols)
        piv.resize(A.num_cols);

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    ValueType *a    = thrust::raw_pointer_cast(&A(0,0));
    IndexType *ipiv = thrust::raw_pointer_cast(&piv[0]);
    lapack_int info = detail::sytrf(order, uplo, n, a, lda, ipiv);

    if( info != 0 )
        throw cusp::runtime_exception("sytrf failed");
}

template<typename Array2d, typename Array1d>
void getrs( const Array2d& A, const Array1d& piv, Array2d& B, char trans )
{
    typedef typename Array2d::value_type ValueType;
    typedef typename Array1d::value_type IndexType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    const ValueType *a    = thrust::raw_pointer_cast(&A(0,0));
    const IndexType *ipiv = thrust::raw_pointer_cast(&piv[0]);

    lapack_int nrhs = B.num_cols;
    lapack_int ldb  = B.pitch;
    ValueType *b    = thrust::raw_pointer_cast(&B(0,0));

    lapack_int info = detail::getrs(order, trans, n, nrhs, a, lda, ipiv, b, ldb);

    if( info != 0 )
        throw cusp::runtime_exception("getrs failed");
}

template<typename Array2d>
void potrs( const Array2d& A, Array2d& B, char uplo )
{
    typedef typename Array2d::value_type ValueType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    const ValueType *a    = thrust::raw_pointer_cast(&A(0,0));

    lapack_int nrhs = B.num_cols;
    lapack_int ldb  = B.pitch;
    ValueType *b    = thrust::raw_pointer_cast(&B(0,0));

    lapack_int info = detail::potrs(order, uplo, n, nrhs, a, lda, b, ldb);

    if( info != 0 )
        throw cusp::runtime_exception("potrs failed");
}

template<typename Array2d, typename Array1d>
void sytrs( const Array2d& A, const Array1d& piv, Array2d& B, char uplo )
{
    typedef typename Array2d::value_type ValueType;
    typedef typename Array1d::value_type IndexType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    const ValueType *a    = thrust::raw_pointer_cast(&A(0,0));
    const IndexType *ipiv = thrust::raw_pointer_cast(&piv[0]);

    lapack_int nrhs = B.num_cols;
    lapack_int ldb  = B.pitch;
    ValueType *b    = thrust::raw_pointer_cast(&B(0,0));

    lapack_int info = detail::sytrs(order, uplo, n, nrhs, a, lda, ipiv, b, ldb);

    if( info != 0 )
        throw cusp::runtime_exception("sytrs failed");
}

template<typename Array2d>
void trtrs( const Array2d& A, Array2d& B, char uplo, char trans, char diag )
{
    typedef typename Array2d::value_type ValueType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    const ValueType *a    = thrust::raw_pointer_cast(&A(0,0));

    lapack_int nrhs = B.num_cols;
    lapack_int ldb  = B.pitch;
    ValueType *b    = thrust::raw_pointer_cast(&B(0,0));

    lapack_int info = detail::trtrs(order, uplo, trans, diag, n, nrhs, a, lda, b, ldb);

    if( info != 0 )
        throw cusp::runtime_exception("trtrs failed");
}

template<typename Array2d>
void trtri( Array2d& A, char uplo, char diag )
{
    typedef typename Array2d::value_type ValueType;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    lapack_int n = A.num_rows;
    lapack_int lda = A.pitch;
    ValueType *a = thrust::raw_pointer_cast(&A(0,0));
    lapack_int info = detail::trtri(order, uplo, diag, n, a, lda);

    if( info != 0 )
        throw cusp::runtime_exception("trtri failed");
}

template<typename Array2d, typename Array1d>
void syev( const Array2d& A, Array1d& eigvals, Array2d& eigvecs, char uplo )
{
    typedef typename Array2d::value_type ValueType;

    if((eigvecs.num_rows != A.num_rows) || (eigvecs.num_cols != A.num_cols))
        eigvals.resize(A.num_rows, A.num_cols);

    if(eigvals.size() != A.num_cols)
        eigvals.resize(A.num_cols);

    cusp::copy(A, eigvecs);

    lapack_int order = Orientation<typename Array2d::orientation>::type;
    char job  = 'V';

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    ValueType *a    = thrust::raw_pointer_cast(&eigvecs(0,0));
    ValueType *w    = thrust::raw_pointer_cast(&eigvals[0]);
    lapack_int info = detail::syev(order, job, uplo, n, a, lda, w);

    if( info != 0 )
        throw cusp::runtime_exception("syev failed");
}

template<typename Array1d1, typename Array1d2, typename Array1d3, typename Array2d>
void stev( const Array1d1& alphas, const Array1d2& betas, Array1d3& eigvals, Array2d& eigvecs, char job )
{
    typedef typename Array2d::value_type ValueType;

    cusp::array1d<ValueType,cusp::host_memory> temp(betas);
    eigvals = alphas;

    lapack_int order = Orientation<typename Array2d::orientation>::type;
    lapack_int n     = alphas.size();
    lapack_int ldz   = n;
    ValueType *a = thrust::raw_pointer_cast(&eigvals[0]);
    ValueType *b = thrust::raw_pointer_cast(&temp[0]);
    ValueType *z = thrust::raw_pointer_cast(&eigvecs(0,0));

    lapack_int info = detail::stev(order, job, n, a, b, z, ldz);

    if( info != 0 )
    {
        printf("stev failure code : %d\n", info);
        throw cusp::runtime_exception("stev failed");
    }
}

template<typename Array1d1, typename Array1d2, typename Array1d3>
void stev( const Array1d1& alphas, const Array1d2& betas, Array1d3& eigvals )
{
    typedef typename Array1d1::value_type ValueType;

    char job  = EvalsOrEvecs<evals>::type;
    cusp::array2d<ValueType,cusp::host_memory> eigvecs;

    stev(alphas, betas, eigvals, eigvecs, job);
}

template<typename Array2d1, typename Array2d2, typename Array1d, typename Array2d3>
void sygv( const Array2d1& A, const Array2d2& B, Array1d& eigvals, Array2d3& eigvecs )
{
    typedef typename Array2d1::value_type ValueType;
    typedef typename Array2d1::orientation Array2dOrientation;

    eigvecs = A;
    cusp::array2d<ValueType,cusp::host_memory,Array2dOrientation> temp(B);

    lapack_int order = Orientation<Array2dOrientation>::type;
    char itype = GenEigOp<gen_op1>::type;
    char job   = EvalsOrEvecs<evecs>::type;
    char uplo  = UpperOrLower<upper>::type;

    lapack_int n    = A.num_rows;
    lapack_int lda  = A.pitch;
    lapack_int ldb  = B.pitch;
    const ValueType *a    = thrust::raw_pointer_cast(&eigvecs(0,0));
    ValueType *b    = thrust::raw_pointer_cast(&temp(0,0));
    ValueType *w    = thrust::raw_pointer_cast(&eigvals[0]);
    lapack_int info = detail::sygv(order, itype, job, uplo, n, a, lda, b, ldb, w);

    if( info != 0 )
    {
        printf("sygv failure code : %d\n", info);
        throw cusp::runtime_exception("sygv failed");
    }
}

template<typename Array2d, typename Array1d>
void gesv( const Array2d& A, Array2d& B, Array1d& pivots )
{
    typedef typename Array1d::value_type IndexType;
    typedef typename Array2d::value_type ValueType;

    Array2d C = A;

    lapack_int order = Orientation<typename Array2d::orientation>::type;

    if(pivots.size() != A.num_cols)
        pivots.resize(A.num_cols);

    lapack_int n    = C.num_rows;
    lapack_int nrhs = B.num_cols;
    lapack_int ldc  = C.pitch;
    lapack_int ldb  = B.pitch;
    ValueType *c    = thrust::raw_pointer_cast(&C(0,0));
    ValueType *b    = thrust::raw_pointer_cast(&B(0,0));
    IndexType *ipiv = thrust::raw_pointer_cast(&pivots[0]);
    lapack_int info = detail::gesv(order, n, nrhs, c, ldc, ipiv, b, ldb);

    if( info != 0 )
        throw cusp::runtime_exception("gesv failed");
}

} // end namespace lapack
} // end namespace cusp

