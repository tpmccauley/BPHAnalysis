#!/bin/sh

export FILE=$1
export TYPE=$2
export FLAG=$3

export DCOM=`find /bin /usr/bin | grep  dirname | head -1`
export BCOM=`find /bin /usr/bin | grep basename | head -1`
export DIR=`${DCOM} $0`
export TOP=`${DCOM} ${DIR}`

export FTMP=/tmp/`${BCOM} ${FILE}``date +%s`

awk -v TYPE=${TYPE} -v FLAG=${FLAG} '((index($0,TYPE)==0)&&(index($0,FLAG)==0)){print $0}' ${FILE} > ${FTMP}

cp ${FTMP} ${FILE}
rm -f ${FTMP}

