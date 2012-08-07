#!/bin/bash

rm -f Makefile
f=Makefile

case "$PHAML_MUMPS" in
   no)
      if [ -e test24.f90 ]
      then
         mv test24.f90 mtest24.f90
      fi
      if [ -e test25.f90 ]
      then
         mv test25.f90 mtest25.f90
      fi ;;
   yes)
      if [ -e mtest24.f90 ]
      then
         mv mtest24.f90 test24.f90
      fi
      if [ -e mtest25.f90 ]
      then
         mv mtest25.f90 test25.f90
      fi ;;
esac

case "$PHAML_SUPERLU" in
   no)
      if [ -e test26.f90 ]
      then
         mv test26.f90 mtest26.f90
      fi ;;
   yes)
      if [ -e mtest26.f90 ]
      then
         mv mtest26.f90 test26.f90
      fi ;;
esac

source ../mkinc/mkhead.inc

for PROG in `ls test*.f90` ;
do
   TESTN=`echo $PROG | sed -e s/\.f90//`
   echo "	$TESTN.exe \\" 1>>$f ;
done
echo "	phaml_slave" 1>>$f

source ../mkinc/mkslave.inc

for PROG in `ls test*.f90` ;
do
   TESTN=`echo $PROG | sed -e s/\.f90//`
   export TESTN
   source ../mkinc/mkmain.inc ;
done

source ../mkinc/mktail.inc
