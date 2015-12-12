#!/bin/sh

export DIR=`dirname $0`

${DIR}/removeFunction.sh getFromBT
${DIR}/removeFunction.sh getFromPC
${DIR}/changeTypedef.sh  reco::PFCandidate
${DIR}/removeInclude.sh "DataFormats/PatCandidates/interface/PackedCandidate.h"

