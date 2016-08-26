#!/bin/sh

export FILE=$1
export TYPE=$2
export FLAG=$3

export DCOM=`find /bin /usr/bin | grep  dirname | head -1`
export BCOM=`find /bin /usr/bin | grep basename | head -1`
export DIR=`${DCOM} $0`
export TOP=`${DCOM} ${DIR}`

export FTMP=/tmp/`${BCOM} ${FILE}``date +%s`

export FOUND=`awk -v TYPE=${TYPE} -v FLAG=${FLAG} 'BEGIN{ff=0;} ((index($0,TYPE)!=0)&&(index($0,FLAG)!=0)){ff=1;} END{print ff;}' ${FILE}`

if [ ${FOUND} -eq "0" ] ; then
  awk -v TYPE=${TYPE} -v FLAG=${FLAG} 'BEGIN {ff=0;} ((ff==0)&&(index($0,"use")!=0)){ff=1;} ((ff==1)&&(index($0,"export")!=0)){ff=2;} ((ff==1)&&(index($0,"library")!=0)){ff=2;} (ff==2){print "<flags "TYPE"=\042"FLAG"\042/>"; ff=3;} {print $0}' ${FILE} > ${FTMP}
  cp ${FTMP} ${FILE}
  rm -f ${FTMP}
fi

