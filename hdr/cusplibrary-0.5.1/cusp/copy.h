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

/*! \file copy.h
 *  \brief Performs (deep) copy operations between containers and views.
 */

#pragma once

#include <cusp/detail/config.h>

#include <thrust/execution_policy.h>

namespace cusp
{

/*! \addtogroup algorithms Algorithms
 *  \addtogroup matrix_algorithms Matrix Algorithms
 *  \ingroup algorithms
 *  \{
 */

/*! \cond */
template <typename DerivedPolicy, typename SourceType, typename DestinationType>
void copy(const thrust::detail::execution_policy_base<DerivedPolicy> &exec,
          const SourceType& src, DestinationType& dst);
/*! \endcond */

/**
 * \brief Copy one array or matrix to another
 *
 * \tparam SourceType Type of the input matrix to copy
 * \tparam DestinationType Type of the output matrix
 *
 * \param src Input matrix to copy
 * \param dst Output matrix created by copying src to dst
 *
 * \note SourceType and DestinationType must have the same format
 * \note DestinationType will be resized as necessary
 *
 * \par Example
 * \code
 * #include <cusp/array1d.h>
 * #include <cusp/print.h>
 *
 * // include cusp copy header file
 * #include <cusp/copy.h>
 *
 * int main()
 * {
 *   // Allocate a array of size 10
 *   cusp::array1d<int,cusp::host_memory> a(10);
 *
 *   // Create a random array with 10 entries
 *   cusp::random_array<int> rand(10);
 *
 *   // Copy random values from rand into a
 *   cusp::copy(rand, a);
 *
 *   // print the contents of a
 *   cusp::print(a);
 * }
 * \endcode
 *
 * \see \p convert
 */
template <typename SourceType, typename DestinationType>
void copy(const SourceType& src, DestinationType& dst);
/*! \}
 */

} // end namespace cusp

#include <cusp/detail/copy.inl>

