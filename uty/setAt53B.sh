#!/bin/sh

export DIR=`dirname $0`
export TOP=`dirname ${DIR}`

eval `find ${TOP} | grep '_53B$' | awk -F_ '{printf("cp "$0" "$1); for(i=2;i<NF;++i) printf("_"$i); print ";"}'`

