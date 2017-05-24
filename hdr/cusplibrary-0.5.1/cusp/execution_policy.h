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

/*! \file cusp/execution_policy.h
 *  \brief Cusp execution policies.
 */

#pragma once

#include <cusp/detail/config.h>

#include <thrust/execution_policy.h>

// #include the host system's execution_policy header
#define __CUSP_HOST_SYSTEM_EXECUTION_POLICY_HEADER <__CUSP_HOST_SYSTEM_ROOT/execution_policy.h>
#include __CUSP_HOST_SYSTEM_EXECUTION_POLICY_HEADER
#undef __CUSP_HOST_SYSTEM_EXECUTION_POLICY_HEADER

// #include the device system's execution_policy.h header
#define __CUSP_DEVICE_SYSTEM_EXECUTION_POLICY_HEADER <__CUSP_DEVICE_SYSTEM_ROOT/execution_policy.h>
#include __CUSP_DEVICE_SYSTEM_EXECUTION_POLICY_HEADER
#undef __CUSP_DEVICE_SYSTEM_EXECUTION_POLICY_HEADER

namespace cusp
{

template<typename DerivedPolicy>
  struct execution_policy
    : thrust::execution_policy<DerivedPolicy>
{};

} // end namespace cusp
