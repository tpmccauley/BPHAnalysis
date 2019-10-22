#!/bin/sh

export DCOM=`find /bin /usr/bin | grep  dirname | head -1`
export BCOM=`find /bin /usr/bin | grep basename | head -1`
export DIR=`${DCOM} $0`
export TOP=`${DCOM} ${DIR}`

export FILE=${TOP}/RecoDecay/interface/BPHTrackReference.h
export FTMP=/tmp/`${BCOM} ${FILE}``date +%s`

export FOUND=`awk 'BEGIN{ff=0;} {cl=index($0,"//"); if(cl==0)cl=length($0);} ((index($0,"dynamic_cast")!=0)&&(index($0,"dynamic_cast")<index($0,"PackedCandidate"))&&(index($0,"PackedCandidate")<cl)){f=1;} ((index($0,"pseudoTrack")!=0)&&(index($0,"pseudoTrack")<cl)){ff=1;} END{print ff;}' ${FILE}`

if [ ${FOUND} -eq "0" ] ; then
  awk 'BEGIN{ff=0;} ((index($0,"static const")!=0)&&(index($0,"getFromPC")!=0)){ff=1;} ((index($0,"return")!=0)&&(ff==1)){print "    const pat::PackedCandidate* pp ="; print "        dynamic_cast<const pat::PackedCandidate*>( &rc );"; print "    if ( pp == 0 ) return 0;"; print "    try {"; print "      const reco::Track* trk = &pp->pseudoTrack();"; print "      return trk;"; print "    }"; print "    catch ( edm::Exception& e ) {"; print "    }"; ff=2;} {print $0;}' ${FILE} > ${FTMP}
  cp ${FTMP} ${FILE}
  rm -f ${FTMP}
fi

