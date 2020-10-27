#!/bin/bash

set -ex

shopt -s expand_aliases
set +u && source ${CMS_PATH}/cmsset_default.sh; set -u
cmsrel ${CMSSW_RELEASE}
cd ${CMSSW_RELEASE}/src
cmsenv
git cms-init --upstream-only

mv ${CI_PROJECT_DIR}/BPHAnalysis ${CMSSW_BASE}/src
scram b
