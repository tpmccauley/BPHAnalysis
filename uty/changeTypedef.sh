#!/bin/sh

export TYPE=$1

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

export FILE=${TOP}/RecoDecay/interface/BPHTrackReference.h
export FTMP=/tmp/`basename ${FILE}``date +%s`

awk -v TYPE=${TYPE} '{line=$0;} (($1=="typedef")&&(substr($3,1,9)=="candidate")){sub($2,TYPE,line);} {print line;}' ${FILE} > ${FTMP}
cp ${FTMP} ${FILE}
rm -f ${FTMP}

