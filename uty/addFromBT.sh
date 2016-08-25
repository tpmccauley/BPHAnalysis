#!/bin/sh

export DCOM=`find /bin /usr/bin | grep  dirname | head -1`
export BCOM=`find /bin /usr/bin | grep basename | head -1`
export DIR=`${DCOM} $0`
export TOP=`${DCOM} ${DIR}`

export FILE=${TOP}/RecoDecay/interface/BPHTrackReference.h
export FTMP=/tmp/`${BCOM} ${FILE}``date +%s`

export FOUND=`awk 'BEGIN{ff=0;} {cl=index($0,"//"); if(cl==0)cl=length($0);} ((index($0,"const reco::Track")!=0)&&(index($0,"const reco::Track")<index($0,"bestTrack"))&&(index($0,"bestTrack")<cl)){ff=1;} END{print ff;}' ${FILE}`

if [ ${FOUND} -eq "0" ] ; then
  awk 'BEGIN{ff=0;} ((index($0,"static const")!=0)&&(index($0,"getFromBT")!=0)){ff=1;} ((index($0,"return")!=0)&&(ff==1)){print "    try {"; print "      const reco::Track* trk = rc.bestTrack();"; print "      return trk;"; print "    }"; print "    catch ( edm::Exception e ) {"; print "    }"; ff=0;} {print $0;}' ${FILE} > ${FTMP}
  cp ${FTMP} ${FILE}
  rm -f ${FTMP}
fi

