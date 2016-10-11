#!/bin/sh

export DSRC=$1
export DEST=$2

cd `dirname $0`
cd ..
export DIR=`pwd`

find ${DIR}/*Decay | egrep '(h|cc|xml|py)$' | grep -v __init__ | awk -F/ -v DSRC=${DSRC} -v DEST=${DEST} '{printf("mkdir -p "); for(i=1;i<NF;++i){s=$i;if(s==DSRC)s=DEST;printf(s);if(i<NF-1)printf("/")} print" ;"; printf("sed s/"DSRC"/"DEST"/ "); for(i=1;i<NF;++i)printf($i"/"); printf($NF" > "); for(i=1;i<=NF;++i){s=$i;if(s==DSRC)s=DEST;printf(s);if(i<NF)printf("/")} print " ;"}'


eval `find ${DIR}/*Decay | egrep '(h|cc|xml|py)$' | grep -v __init__ | awk -F/ -v DSRC=${DSRC} -v DEST=${DEST} '{printf("mkdir -p "); for(i=1;i<NF;++i){s=$i;if(s==DSRC)s=DEST;printf(s);if(i<NF-1)printf("/")} print" ;"; printf("sed s/"DSRC"/"DEST"/ "); for(i=1;i<NF;++i)printf($i"/"); printf($NF" > "); for(i=1;i<=NF;++i){s=$i;if(s==DSRC)s=DEST;printf(s);if(i<NF)printf("/")} print " ;"}'`


