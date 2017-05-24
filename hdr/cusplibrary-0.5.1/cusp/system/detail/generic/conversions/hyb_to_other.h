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


#pragma once

#include <cusp/sort.h>

#include <cusp/system/detail/generic/conversions/ell_to_other.h>

#include <thrust/copy.h>

namespace cusp
{
namespace system
{
namespace detail
{
namespace generic
{

template <typename DerivedPolicy, typename SourceType, typename DestinationType>
typename enable_if_same_system<SourceType,DestinationType>::type
convert(thrust::execution_policy<DerivedPolicy>& exec,
        const SourceType& src,
        DestinationType& dst,
        cusp::hyb_format&,
        cusp::coo_format&)
{
    typedef typename SourceType::const_coo_view_type  CooViewType;

    CooViewType src_coo(src);
    cusp::copy(src_coo, dst);
}

} // end namespace generic
} // end namespace detail
} // end namespace system
} // end namespace cusp
