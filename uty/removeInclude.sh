#!/bin/sh

export FILE=$1
export INCF=$2
export INCN=`echo ${INCF} | awk -F/ '{print $NF}' | awk -F. '{print $1}'`

export BCOM=`find /bin /usr/bin | grep basename | head -1`
export FTMP=/tmp/`${BCOM} ${FILE}``date +%s`

awk -v INCN="${INCN}" '{cl=index($0,"//"); if(cl==0)cl=length($0);} ((index($0,"#include")==0)||(index($0,INCN)==0)||(index($0,INCN)>cl)){print $0;}' ${FILE} > ${FTMP}
cp ${FTMP} ${FILE}
rm -f ${FTMP}

