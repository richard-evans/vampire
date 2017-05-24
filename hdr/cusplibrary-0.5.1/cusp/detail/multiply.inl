/*
 *  Copyright 2008-2014 NVIDIA Corporation
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

#include <thrust/detail/config.h>
#include <thrust/system/detail/generic/select_system.h>

#include <cusp/multiply.h>

#include <cusp/detail/type_traits.h>
#include <cusp/system/detail/adl/multiply.h>
#include <cusp/system/detail/generic/multiply.h>

namespace cusp
{

template <typename DerivedPolicy,
         typename LinearOperator,
         typename MatrixOrVector1,
         typename MatrixOrVector2>
void multiply(const thrust::detail::execution_policy_base<DerivedPolicy> &exec,
              const LinearOperator&  A,
              const MatrixOrVector1& B,
              MatrixOrVector2& C)
{
    using cusp::system::detail::generic::multiply;

    typedef typename LinearOperator::format  Format1;
    typedef typename MatrixOrVector1::format Format2;
    typedef typename MatrixOrVector2::format Format3;

    Format1 format1;
    Format2 format2;
    Format3 format3;

    multiply(thrust::detail::derived_cast(thrust::detail::strip_const(exec)),
             A, B, C,
             format1, format2, format3);
}

template <typename LinearOperator,
         typename MatrixOrVector1,
         typename MatrixOrVector2>
void multiply(const LinearOperator&  A,
              const MatrixOrVector1& B,
              MatrixOrVector2& C)
{
    using thrust::system::detail::generic::select_system;

    typedef typename LinearOperator::memory_space  System1;
    typedef typename MatrixOrVector1::memory_space System2;
    typedef typename MatrixOrVector2::memory_space System3;

    System1 system1;
    System2 system2;
    System3 system3;

    cusp::multiply(select_system(system1,system2,system3), A, B, C);
}

template <typename DerivedPolicy,
         typename LinearOperator,
         typename MatrixOrVector1,
         typename MatrixOrVector2,
         typename UnaryFunction,
         typename BinaryFunction1,
         typename BinaryFunction2>
typename thrust::detail::disable_if_convertible<UnaryFunction,cusp::known_format,void>::type
multiply(const thrust::detail::execution_policy_base<DerivedPolicy> &exec,
         const LinearOperator&  A,
         const MatrixOrVector1& B,
         MatrixOrVector2& C,
         UnaryFunction  initialize,
         BinaryFunction1 combine,
         BinaryFunction2 reduce)
{
    using cusp::system::detail::generic::multiply;

    typedef typename LinearOperator::format  Format1;
    typedef typename MatrixOrVector1::format Format2;
    typedef typename MatrixOrVector2::format Format3;

    Format1 format1;
    Format2 format2;
    Format3 format3;

    multiply(thrust::detail::derived_cast(thrust::detail::strip_const(exec)),
             A, B, C,
             initialize, combine, reduce,
             format1, format2, format3);
}

template <typename LinearOperator,
         typename MatrixOrVector1,
         typename MatrixOrVector2,
         typename UnaryFunction,
         typename BinaryFunction1,
         typename BinaryFunction2>
void multiply(const LinearOperator&  A,
              const MatrixOrVector1& B,
              MatrixOrVector2& C,
              UnaryFunction  initialize,
              BinaryFunction1 combine,
              BinaryFunction2 reduce)
{
    using thrust::system::detail::generic::select_system;

    typedef typename LinearOperator::memory_space  System1;
    typedef typename MatrixOrVector1::memory_space System2;
    typedef typename MatrixOrVector2::memory_space System3;

    System1 system1;
    System2 system2;
    System3 system3;

    cusp::multiply(select_system(system1,system2,system3), A, B, C,
                   initialize, combine, reduce);
}

template <typename DerivedPolicy,
         typename LinearOperator,
         typename Vector1,
         typename Vector2,
         typename Vector3,
         typename BinaryFunction1,
         typename BinaryFunction2>
void generalized_spmv(const thrust::detail::execution_policy_base<DerivedPolicy> &exec,
                      const LinearOperator&  A,
                      const Vector1& x,
                      const Vector2& y,
                      Vector3& z,
                      BinaryFunction1 combine,
                      BinaryFunction2 reduce)
{
    using cusp::system::detail::generic::generalized_spmv;

    typedef typename LinearOperator::format  Format1;
    typedef typename Vector1::format         Format2;
    typedef typename Vector2::format         Format3;
    typedef typename Vector3::format         Format4;

    Format1 format1;
    Format2 format2;
    Format3 format3;
    Format4 format4;

    generalized_spmv(thrust::detail::derived_cast(thrust::detail::strip_const(exec)),
                     A, x, y, z,
                     combine, reduce,
                     format1, format2, format3, format4);
}

template <typename LinearOperator,
         typename Vector1,
         typename Vector2,
         typename Vector3,
         typename BinaryFunction1,
         typename BinaryFunction2>
void generalized_spmv(const LinearOperator&  A,
                      const Vector1& x,
                      const Vector2& y,
                      Vector3& z,
                      BinaryFunction1 combine,
                      BinaryFunction2 reduce)
{
    using thrust::system::detail::generic::select_system;

    typedef typename LinearOperator::memory_space  System1;
    typedef typename Vector1::memory_space System2;
    typedef typename Vector2::memory_space System3;
    typedef typename Vector3::memory_space System4;

    System1 system1;
    System2 system2;
    System3 system3;
    System4 system4;

    cusp::generalized_spmv(select_system(system1,system2,system3,system4), A, x, y, z, combine, reduce);
}

} // end namespace cusp

