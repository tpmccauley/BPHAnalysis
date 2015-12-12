#!/bin/sh

export DIR=`dirname $0`

${DIR}/addInclude.sh    DataFormats/PatCandidates/interface/PackedCandidate.h
${DIR}/changeTypedef.sh pat::PackedCandidate
${DIR}/addFromBT.sh
${DIR}/addFromPC.sh

