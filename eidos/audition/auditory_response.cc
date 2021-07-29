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

#include "eidos/audition/auditory_response.h"

#include <cstdint>

#include "cnpy.h"
#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

// Name of the array containing the center frequencies.
constexpr char kCenterFrequenciesName[] = "CENTER_FREQUENCIES";

}  // namespace

void AuditoryResponse::WriteNpy(const std::string &filename) const {
  GOOGLE_CHECK_LT(0, outputs_.size());
  const auto &response = outputs_.back().second;
  const uint32_t num_rows = response.rows();
  const uint32_t num_cols = response.cols();
  const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      row_mat = response.cast<float>();
  const std::vector<std::size_t> shape = {static_cast<std::size_t>(num_rows),
                                          static_cast<std::size_t>(num_cols)};
  const float *data = row_mat.data();
  cnpy::npy_save(filename, data, shape, "w");
}

void AuditoryResponse::WriteNpz(const std::string &filename) const {
  std::string write_mode = "w";
  GOOGLE_LOG(INFO) << "Writing " << outputs_.size() << " outputs to "
                   << filename;
  for (uint32_t i = 0; i < outputs_.size(); ++i) {
    const Eigen::ArrayXXd &response = outputs_[i].second;
    const std::size_t num_rows = response.rows();
    const std::size_t num_cols = response.cols();
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        row_mat = response.cast<float>();
    const std::vector<std::size_t> shape = { num_rows, num_cols };
    const float *data = row_mat.data();
    const std::string output_name = AuditoryOutputType_Name(outputs_[i].first);
    GOOGLE_LOG(INFO) << "Writing " << output_name;
    cnpy::npz_save(filename, output_name, data, shape, write_mode);
    write_mode = "a";
  }

  // Save center frequencies as a separate vector.
  const auto &center_frequencies = channel_properties_.center_frequencies;
  const std::size_t num_channels = center_frequencies.size();
  GOOGLE_CHECK_LT(0, num_channels);
  cnpy::npz_save(filename, kCenterFrequenciesName, center_frequencies.data(),
                 { num_channels }, write_mode);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
