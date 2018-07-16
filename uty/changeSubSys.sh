#!/bin/sh

cd `dirname $0`
cd ..
export DIR=`pwd`

if [ $# -ge "2" ] ; then
  export DSRC=$1
  export DEST=$2
else
  export DSRC=`echo ${DIR} | awk -F/ '{print $NF}'`
  export DEST=$1
fi

cd ..

export DLIST=`ls -1d ${DSRC}/* | grep -v uty | awk '{for(i=1;i<NF;++i)printf($i" "); print $NF}'`
find ${DLIST} | egrep '(h|cc|C|xml|py)$' | grep -v __init__ | awk -F/ -v DSRC=${DSRC} -v DEST=${DEST} '(NF>1){p=""; for(i=2;i<NF;++i)p=p"/"$i; s=DSRC""p; d=DEST""p; print "mkdir -p "d" ; sed s/"DSRC"/"DEST"/ "s"/"$NF" > "d"/"$NF" ;"}'
eval `find ${DLIST} | egrep '(h|cc|C|xml|py)$' | grep -v __init__ | awk -F/ -v DSRC=${DSRC} -v DEST=${DEST} '(NF>1){p=""; for(i=2;i<NF;++i)p=p"/"$i; s=DSRC""p; d=DEST""p; print "mkdir -p "d" ; sed s/"DSRC"/"DEST"/ "s"/"$NF" > "d"/"$NF" ;"}'`

cp -r --preserve=all ${DSRC}/uty ${DEST}/uty
