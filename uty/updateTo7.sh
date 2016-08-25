#!/bin/sh

export DCOM=`find /bin /usr/bin | grep dirname | head -1`
export DIR=`${DCOM} $0`
export TOP=`${DCOM} ${DIR}`

${DIR}/addInclude.sh    ${TOP}/RecoDecay/interface/BPHTrackReference.h \
                        DataFormats/PatCandidates/interface/PackedCandidate.h
${DIR}/addFromBT.sh
${DIR}/addFromPC.sh
${DIR}/changeTypedef.sh ${TOP}/RecoDecay/interface/BPHTrackReference.h \
                        pat::PackedCandidate candidate
${DIR}/changeTypedef.sh ${TOP}/SpecificDecay/plugins/BPHWriteSpecificDecay.h \
                        std::unique_ptr\<pat::CompositeCandidateCollection\> ccc_pointer

eval `find ${TOP}/RecoDecay     | grep 'BuildFile.xml$' | awk -v DIR=${DIR} '{print DIR"/removeFlags.sh "$0" CXXFLAGS -Wno-delete-non-virtual-dtor ;"}'`
eval `find ${TOP}/SpecificDecay | grep 'BuildFile.xml$' | awk -v DIR=${DIR} '{print DIR"/removeFlags.sh "$0" CXXFLAGS -Wno-delete-non-virtual-dtor ;"}'`

