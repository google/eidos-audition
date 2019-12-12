// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Simple checks that the Boost dependencies work.

#include <utility>

#include "absl/strings/str_join.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/regex.hpp"
#include "gtest/gtest.h"

using namespace boost;

namespace eidos {

TEST(BoostTest, CheckRegex) {
  const char *kRe = "^[[:space:]]*";
  const regex expression(kRe);
  EXPECT_TRUE(regex_match("   ", expression));
  EXPECT_FALSE(regex_match("abc", expression));
}

TEST(BoostTest, GraphTest) {
  // Create a typedef for the Graph type.
  typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

  // Make convenient labels for the vertices.
  enum { A, B, C, D, E, N };
  const int num_vertices = N;

  // Writing out the edges in the graph.
  typedef std::pair<int, int> Edge;
  Edge edge_array[] = {Edge(A, B), Edge(A, D), Edge(C, A), Edge(D, C),
                       Edge(C, E), Edge(B, D), Edge(D, E)};
  const int num_edges = sizeof(edge_array) / sizeof(edge_array[0]);

  // Declare a graph object.
  Graph g(num_vertices);

  // Add the edges to the graph object.
  for (int i = 0; i < num_edges; ++i) {
    add_edge(edge_array[i].first, edge_array[i].second, g);
  }

  // Get the property map for vertex indices.
  typedef property_map<Graph, vertex_index_t>::type IndexMap;
  IndexMap index = get(vertex_index, g);

  typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> vp;
  std::vector<int> vertix_indices;
  vertix_indices.reserve(N);
  for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
    vertix_indices.push_back(index[*vp.first]);
  }
  EXPECT_EQ("0 1 2 3 4", absl::StrJoin(vertix_indices, " "));
}

}  // namespace eidos

// Local Variables:
// mode: c++
// End:
