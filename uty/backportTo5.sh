#!/bin/sh

export DCOM=`find /bin /usr/bin | grep dirname | head -1`
export DIR=`${DCOM} $0`
export TOP=`${DCOM} ${DIR}`

${DIR}/removeFunction.sh getFromBT
${DIR}/removeFunction.sh getFromPC
${DIR}/changeTypedef.sh ${TOP}/RecoDecay/interface/BPHTrackReference.h \
                        reco::PFCandidate candidate
${DIR}/removeInclude.sh ${TOP}/RecoDecay/interface/BPHTrackReference.h \
                        "DataFormats/PatCandidates/interface/PackedCandidate.h"
${DIR}/changeTypedef.sh ${TOP}/SpecificDecay/plugins/BPHWriteSpecificDecay.h \
                        std::auto_ptr\<pat::CompositeCandidateCollection\> ccc_pointer

eval `find ${TOP}/RecoDecay     | grep 'BuildFile.xml$' | awk -v DIR=${DIR} '{print DIR"/addFlags.sh "$0" CXXFLAGS -Wno-delete-non-virtual-dtor ;"}'`
eval `find ${TOP}/SpecificDecay | grep 'BuildFile.xml$' | awk -v DIR=${DIR} '{print DIR"/addFlags.sh "$0" CXXFLAGS -Wno-delete-non-virtual-dtor ;"}'`

rm -f ${TOP}/SpecificDecay/src/classes*

