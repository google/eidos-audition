# Copyright 2019 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Makefile for building DSAM test drivers.

CC = g++
CFLAGS = -I.. -I/usr/local/include/dsam
DEFINES = -DHAVE_CONFIG_H -DDSAM_VERSION_INT=2008044 -DDSAM_USE_UNICODE=1
LFLAGS = -L/usr/local/lib -ldsam -lfftw3 -lsndfile -lm

all: auditory_nerve_test hair_cell_synapse_test

auditory_nerve_test: auditory_nerve_test.cc input_data_file.params bm_gammatone.params ihc_meddis_86.params anf_spikes_binomial.params
	$(CC) $(CFLAGS) $(DEFINES) -O2 auditory_nerve_test.cc -o auditory_nerve_test $(LFLAGS)

hair_cell_synapse_test: hair_cell_synapse_test.cc input_data_file.params bm_gammatone.params ihc_meddis_86.params
	$(CC) $(CFLAGS) $(DEFINES) -O2 hair_cell_synapse_test.cc -o hair_cell_synapse_test $(LFLAGS)

clean:
	rm -f auditory_nerve_test hair_cell_synapse_test output.dat *.diag *~
