#!/bin/sh

export INCF=$1
export INCN=`echo ${INCF} | awk -F/ '{print $NF}' | awk -F. '{print $1}'`

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

export FILE=${TOP}/RecoDecay/interface/BPHTrackReference.h
export FTMP=/tmp/`basename ${FILE}``date +%s`

awk -v INCN="${INCN}" '{cl=index($0,"//"); if(cl==0)cl=length($0);} ((index($0,"#include")==0)||(index($0,INCN)==0)||(index($0,INCN)>cl)){print $0;}' ${FILE} > ${FTMP}
cp ${FTMP} ${FILE}
rm -f ${FTMP}

