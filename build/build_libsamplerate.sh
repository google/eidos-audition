#!/bin/bash
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

# This utility fetches the current snapshot of libresample,
# configures and builds it.

set -o pipefail

echo "=== Building libresample from the current HEAD ..."

if [ $# -ne 1 ] ; then
    echo "Usage: $0 <build_dir>"
    exit 1
fi

BUILD_DIR=$1
ORIGINAL_DIR=$PWD

# Create build directory if it doesn't exist.
TMP_BUILD_DIR="/tmp/build_libresample"
mkdir -p "${TMP_BUILD_DIR}" > /dev/null 2>&1
if [ $? -ne 0 ] ; then
    echo "Failed to create build directory"
    exit 1
fi

cd $TMP_BUILD_DIR

# Fetch the tarball at HEAD.
TARBALL="master.tar.gz"

echo "=== Fetching the source ..."
WGET=$(which wget)
if [ $? -ne 0 ] ; then
    echo "wget not found"
    exit 1
fi
"${WGET}" "https://github.com/erikd/libsamplerate/archive/${TARBALL}"
if [ $? -ne 0 ] ; then
    echo "Failed to fetch the tarball"
    exit 1
fi

# Unpack.
SOURCE_DIR="libsamplerate-master"
echo "=== Unpacking the source ..."
TAR=$(which tar)
if [ $? -ne 0 ] ; then
    echo "tar not found"
    exit 1
fi
"${TAR}" xfz "${TARBALL}"
if [ $? -ne 0 ] ; then
    echo "Failed to unpack archive"
    exit 1
fi

# Configure and build.
#
# Note that the shared library creation is disabled. This is to enforce
# the compatibility of builds between Linux and Darwin platforms.
OUTPUT_DIR="${TMP_BUILD_DIR}/outputs"
cd "${SOURCE_DIR}"
if [ ! -f "autogen.sh" ] ; then
    echo "Autoconf bootstrapper not found"
    exit 1
fi
./autogen.sh
if [ $? -ne 0 ] ; then
    echo "Failed to bootstrap autoconf"
    exit 1
fi
# We disable gmp because for some reason the xcode >10.0 linker fails
# to find it.
./configure --prefix="${OUTPUT_DIR}" --disable-shared
if [ $? -ne 0 ] ; then
    echo "Configure failed"
    exit 1
fi
mkdir -p "${OUTPUT_DIR}"
if [ $? -ne 0 ] ; then
    echo "Failed to create ${OUTPUT_DIR}"
    exit 1
fi
make -j 5
if [ $? -ne 0 ] ; then
    echo "make failed"
    exit 1
fi
make install
if [ $? -ne 0 ] ; then
    echo "make install failed"
    exit 1
fi
echo "=== Build complete. Artifacts can be found in: ${OUTPUT_DIR}"

# Copy the artifacts to the supplied Bazel build dir.
cd "${ORIGINAL_DIR}"
echo "=== Copying artifacts to \"${BUILD_DIR}\""
cp "${OUTPUT_DIR}/lib/libsamplerate.a" "${BUILD_DIR}"
if [ $? -ne 0 ] ; then
    echo "Failed to copy library artifacts to Bazel directory"
    exit 1
fi
cp "${OUTPUT_DIR}/include/"*.h "${BUILD_DIR}"
if [ $? -ne 0 ] ; then
    echo "Failed to copy header artifacts to Bazel directory"
    exit 1
fi

echo "=== Cleaning up temporary build space"
rm -rf "${TMP_BUILD_DIR}"

exit 0
