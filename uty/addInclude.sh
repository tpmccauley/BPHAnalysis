#!/bin/sh

export INCF=$1
export INCN=`echo ${INCF} | awk -F/ '{print $NF}' | awk -F. '{print $1}'`

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

export FILE=${TOP}/RecoDecay/interface/BPHTrackReference.h
export FTMP=/tmp/`basename ${FILE}``date +%s`

export FOUND=`awk -v INCN=${INCN} 'BEGIN{ff=0;} {cl=index($0,"//"); if(cl==0)cl=length($0);} ((index($0,"#include")!=0)&&(index($0,INCN)!=0)&&(index($0,INCN)<cl)){ff=1;} END{print ff;}' ${FILE}`

if [ ${FOUND} -eq "0" ] ; then
  awk -v INCF=${INCF} '((index($0,"#include")!=0)&&(index($0,"PFCandidate")!=0)){print "#include \042"INCF"\042";} {print $0;}' ${FILE} > ${FTMP}
  cp ${FTMP} ${FILE}
  rm -f ${FTMP}
fi

