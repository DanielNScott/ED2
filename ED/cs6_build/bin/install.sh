#!/bin/sh
if [ 'x'${1} == 'xclean' ]
then
  make OPT=cs6${2} clean
else
  make OPT=cs6${2}
fi
