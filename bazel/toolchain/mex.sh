#!/bin/bash
#
# Copyright 2015 The Bazel Authors. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# OS X relpath is not really working. This is a wrapper script around gcc
# to simulate relpath behavior.
#
# This wrapper uses install_name_tool to replace all paths in the binary
# (bazel-out/.../path/to/original/library.so) by the paths relative to
# the binary. It parses the command line to behave as rpath is supposed
# to work.
#
# See https://blogs.oracle.com/dipol/entry/dynamic_libraries_rpath_and_mac
# on how to set those paths for Mach-O binaries.
#
set -eu

GCC=/Applications/MATLAB_R2022a.app/bin/mex
OUTPUT=
DEPENDENCY_FILE=

declare -a commands

function parse_option() {
    local -r opt="$1"
    echo Found out "$opt"
    if [[ "${OUTPUT}" = "1" ]]; then
        OUTPUT=$(dirname "$opt")
        commands+=("$OUTPUT")
    elif [[ "${DEPENDENCY_FILE}" = "1" ]]; then
        DEPENDENCY_FILE="$opt"
    elif [[ "$opt" = "-o" ]]; then
        # output is coming
        OUTPUT=1
        commands+=("-outdir")
    elif [[ "$opt" = "-dependency_file" ]]; then
        # output is coming
        DEPENDENCY_FILE=1
    else
        commands+=("$opt")
    fi
}

# let parse the option list
for i in "$@"; do
    if [[ "$i" = @* ]]; then
        while IFS= read -r opt
        do
            parse_option "$opt"
        done < "${i:1}" || exit 1
    else
        parse_option "$i"
    fi
done

# Call gcc
echo "Calling MATLAB_R2022a"
echo "${GCC}" "${commands[@]}"

set -v
${GCC} "${commands[@]}"

PWD=$(pwd)
echo sed -i '.bak' "s#${PWD}/##" ${DEPENDENCY_FILE}
sed -i '.bak' "s#${PWD}/##" ${DEPENDENCY_FILE}
# TODO correct dependency file by stripping PWD