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

// Simple DSAM-based test driver for hair cell models.
//
// The parameter files are copied from either DSAM examples or SAMS tutorials.

#include <DSAM.h>
#include <iostream>

int main() {
  // Initialize modules.
  EarObjectPtr stimulus = Init_EarObject(wxT("DataFile_In"));
  EarObjectPtr bm_filter = Init_EarObject(wxT("BM_gammaT"));
  EarObjectPtr hair_cell = Init_EarObject(wxT("IHC_Meddis86"));
  EarObjectPtr anf_spikes = Init_EarObject(wxT("An_SG_Binomial"));

  // Setup connections between modules.
  ConnectOutSignalToIn_EarObject(stimulus, bm_filter);
  ConnectOutSignalToIn_EarObject(bm_filter, hair_cell);
  ConnectOutSignalToIn_EarObject(hair_cell, anf_spikes);

  // Read module configurations.
  ReadPars_ModuleMgr(stimulus, wxT("input_data_file.params"));
  PrintPars_ModuleMgr(stimulus);
  ReadPars_ModuleMgr(bm_filter, wxT("bm_gammatone.params"));
  PrintPars_ModuleMgr(bm_filter);
  ReadPars_ModuleMgr(hair_cell, wxT("ihc_meddis_86.params"));
  PrintPars_ModuleMgr(hair_cell);
  ReadPars_ModuleMgr(anf_spikes, wxT("anf_spikes_binomial.params"));
  PrintPars_ModuleMgr(anf_spikes);

  // Run the simulation.
  RunProcess_ModuleMgr(stimulus);
  RunProcess_ModuleMgr(bm_filter);
  RunProcess_ModuleMgr(hair_cell);
  RunProcess_ModuleMgr(anf_spikes);

  WriteOutSignal_DataFile(wxT("output.dat"), anf_spikes);
  FreeAll_EarObject();
  return 0;
}

// Local Variables:
// mode: c++
// End:
