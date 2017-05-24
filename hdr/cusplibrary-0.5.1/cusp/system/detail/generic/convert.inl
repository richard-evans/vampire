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


#include <cusp/copy.h>
#include <cusp/detail/format.h>

#include <cusp/format_utils.h>

#include <cusp/system/detail/generic/conversions/array_to_other.h>
#include <cusp/system/detail/generic/conversions/coo_to_other.h>
#include <cusp/system/detail/generic/conversions/csr_to_other.h>
#include <cusp/system/detail/generic/conversions/dia_to_other.h>
#include <cusp/system/detail/generic/conversions/ell_to_other.h>
#include <cusp/system/detail/generic/conversions/hyb_to_other.h>
#include <cusp/system/detail/generic/conversions/permutation_to_other.h>

namespace cusp
{
namespace system
{
namespace detail
{
namespace generic
{
namespace detail
{

template <typename DerivedPolicy,
         typename SourceType,
         typename DestinationType,
         typename Format>
void convert(thrust::execution_policy<DerivedPolicy>& exec,
             const SourceType& src,
             DestinationType& dst,
             Format&,
             Format&)
{
    cusp::copy(src, dst);
}

template <typename DerivedPolicy,
         typename SourceType,
         typename DestinationType,
         typename Format1,
         typename Format2>
typename enable_if_same_system<SourceType,DestinationType>::type
convert(thrust::execution_policy<DerivedPolicy>& exec,
        const SourceType& src,
        DestinationType& dst,
        Format1&,
        Format2&)
{
    // convert src -> coo_matrix -> dst
    typedef typename SourceType::container ContainerType;
    typename cusp::detail::as_coo_type<ContainerType>::type tmp;

    cusp::convert(exec, src, tmp);
    cusp::convert(exec, tmp, dst);
}

template <typename DerivedPolicy,
         typename SourceType,
         typename DestinationType,
         typename Format1,
         typename Format2>
typename enable_if_different_system<SourceType,DestinationType>::type
convert(thrust::execution_policy<DerivedPolicy>& exec,
        const SourceType& src,
        DestinationType& dst,
        Format1&,
        Format2&)
{
    typedef typename SourceType::memory_space MemorySpace;
    typedef typename DestinationType::format  DestFormat;
    typedef typename cusp::detail::as_matrix_type<SourceType,MemorySpace,DestFormat>::type SrcDestType;

    SrcDestType tmp;

    cusp::convert(src, tmp);
    cusp::copy(tmp, dst);
}

} // end namespace detail

template <typename DerivedPolicy,
         typename SourceType,
         typename DestinationType,
         typename Format1,
         typename Format2>
void convert(thrust::execution_policy<DerivedPolicy>& exec,
             const SourceType& src,
             DestinationType& dst,
             Format1& format1,
             Format2& format2)
{
    detail::convert(exec, src, dst, format1, format2);
}

} // end namespace generic
} // end namespace detail
} // end namespace system
} // end namespace cusp
