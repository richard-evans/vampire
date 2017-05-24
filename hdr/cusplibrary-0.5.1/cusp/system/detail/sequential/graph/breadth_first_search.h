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

#include <deque>

#include <cusp/detail/config.h>
#include <cusp/detail/format.h>

#include <cusp/array1d.h>
#include <cusp/exception.h>

#include <cusp/system/detail/sequential/execution_policy.h>

namespace cusp
{
namespace system
{
namespace detail
{
namespace sequential
{

template<typename DerivedPolicy, typename MatrixType, typename ArrayType>
void breadth_first_search(sequential::execution_policy<DerivedPolicy>& exec,
                          const MatrixType& G,
                          const typename MatrixType::index_type src,
                          ArrayType& labels,
                          const bool mark_levels,
                          csr_format)
{
    typedef typename MatrixType::index_type VertexId;

    ArrayType predecessors;

    // initialize predecessor array
    if(!mark_levels) predecessors.resize(G.num_rows);

    // initialize distances
    for (size_t i = 0; i < G.num_rows; i++)
    {
        labels[i] = -1;
    }

    labels[src] = 0;
    VertexId search_depth = 0;

    // Initialize queue for managing previously-discovered nodes
    std::deque<VertexId> frontier;
    frontier.push_back(src);

    //
    // Perform BFS
    //
    while (!frontier.empty()) {

        // Dequeue node from frontier
        VertexId dequeued_node = frontier.front();
        frontier.pop_front();
        VertexId neighbor_dist = labels[dequeued_node] + 1;

        // Locate adjacency list
        int edges_begin = G.row_offsets[dequeued_node];
        int edges_end = G.row_offsets[dequeued_node + 1];

        for (int edge = edges_begin; edge < edges_end; edge++) {

            // Lookup neighbor and enqueue if undiscovered
            VertexId neighbor = G.column_indices[edge];
            if (neighbor == -1) continue;

            if (labels[neighbor] == -1)
            {
                labels[neighbor] = neighbor_dist;

                if (search_depth < neighbor_dist)
                {
                    search_depth = neighbor_dist;
                }

                if(!mark_levels) predecessors[neighbor] = dequeued_node;
                frontier.push_back(neighbor);
            }
        }
    }
    search_depth++;

    // if predecessors are needed then copy into outgoing array
    if(!mark_levels) labels = predecessors;
}

} // end namespace sequential
} // end namespace detail
} // end namespace system

// hack until ADL is operational
using cusp::system::detail::sequential::breadth_first_search;

} // end namespace cusp
