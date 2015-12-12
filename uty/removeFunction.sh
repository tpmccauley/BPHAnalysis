#!/bin/sh

export FUNC=$1

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

export FILE=${TOP}/RecoDecay/interface/BPHTrackReference.h
export FTMP=/tmp/`basename ${FILE}``date +%s`

awk -v FUNC=${FUNC} 'BEGIN{ff=0;} ((ff==0)||(index($1,"//")>0)){print $0;} ((index($0,"static const")!=0)&&(index($0,FUNC)!=0)){ff=1;cc=0;} (ff>0){line=$0; while(index(line,"{")>0){line=substr(line,index(line,"{")+1);++cc;}; while(index(line,"}")>0){line=substr(line,index(line,"}")+1);--cc;};} (cc>0){ff=2} ((ff==2)&&(cc==0)){print "    return 0;"; print "  }";ff=0;}' ${FILE} > ${FTMP}
cp ${FTMP} ${FILE}
rm -f ${FTMP}

