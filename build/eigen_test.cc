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

// Simple checks that Eigen behaves.

#include "Eigen/Dense"
#include "gtest/gtest.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace eidos {

TEST(EigenTest, Trivial) {
  const MatrixXd m = MatrixXd::Identity(3, 3);
  EXPECT_EQ(3, m.rows());
  EXPECT_EQ(3, m.cols());
  EXPECT_EQ(1.0, m(0, 0));
  EXPECT_EQ(1.0, m(1, 1));
  EXPECT_EQ(1.0, m(2, 2));
  VectorXd v(3);
  v << 1, 2, 3;
  EXPECT_EQ(1, v(0));
  EXPECT_EQ(2, v(1));
  EXPECT_EQ(3, v(2));
}

TEST(EigenTest, CheckReverse) {
  MatrixXd m(3, 2);
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_EQ(1, m(0, 0));
  EXPECT_EQ(2, m(0, 1));
  EXPECT_EQ(3, m(1, 0));
  EXPECT_EQ(4, m(1, 1));
  EXPECT_EQ(5, m(2, 0));
  EXPECT_EQ(6, m(2, 1));
  m.colwise().reverseInPlace();
  EXPECT_EQ(5, m(0, 0));
  EXPECT_EQ(6, m(0, 1));
  EXPECT_EQ(3, m(1, 0));
  EXPECT_EQ(4, m(1, 1));
  EXPECT_EQ(1, m(2, 0));
  EXPECT_EQ(2, m(2, 1));
}

}  // namespace eidos

// Local Variables:
// mode: c++
// End:
