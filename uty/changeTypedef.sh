#!/bin/sh

export FILE=$1
export TYPE=$2
export NAME=$3

export BCOM=`find /bin /usr/bin | grep basename | head -1`
export FTMP=/tmp/`${BCOM} ${FILE}``date +%s`

awk -v NAME=${NAME} -v TYPE=${TYPE} '{line=$0;} (($1=="typedef")&&(substr($3,1,length(NAME))==NAME)){sub($2,TYPE,line);} {print line;}' ${FILE} > ${FTMP}
cp ${FTMP} ${FILE}
rm -f ${FTMP}

