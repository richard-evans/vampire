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

#include <cusp/detail/config.h>
#include <cusp/detail/format.h>
#include <cusp/system/detail/sequential/execution_policy.h>

#include <cusp/system/detail/sequential/multiply/coo_spmv.h>
#include <cusp/system/detail/sequential/multiply/csr_spmv.h>
#include <cusp/system/detail/sequential/multiply/dia_spmv.h>
#include <cusp/system/detail/sequential/multiply/ell_spmv.h>
#include <cusp/system/detail/sequential/multiply/hyb_spmv.h>

#include <cusp/system/detail/sequential/multiply/array2d_mv.h>
#include <cusp/system/detail/sequential/multiply/array2d_mm.h>

#include <cusp/system/detail/sequential/multiply/csr_spgemm.h>
#include <cusp/system/detail/sequential/multiply/coo_spgemm.h>

namespace cusp
{
namespace system
{
namespace detail
{
namespace sequential
{

} // end namespace sequential
} // end namespace detail
} // end namespace system

// hack until ADL is operational
using cusp::system::detail::sequential::multiply;

} // end namespace cusp

