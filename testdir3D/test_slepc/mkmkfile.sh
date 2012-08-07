#!/bin/bash

rm -f Makefile
f=Makefile

case "$PHAML_ARPACK" in
   no)
      if [ -e test16.f90 ]
      then
         mv test16.f90 mtest16.f90
      fi
      if [ -e test17.f90 ]
      then
         mv test17.f90 mtest17.f90
      fi
      if [ -e test18.f90 ]
      then
         mv test18.f90 mtest18.f90
      fi ;;
   yes)
      if [ -e mtest16.f90 ]
      then
         mv mtest16.f90 test16.f90
      fi
      case "$PHAML_MUMPS" in
         no)
            if [ -e test17.f90 ]
            then
               mv test17.f90 mtest17.f90
            fi ;;
         yes)
            if [ -e mtest17.f90 ]
            then
               mv mtest17.f90 test17.f90
            fi ;;
      esac
      case "$PHAML_SUPERLU" in
         no)
            if [ -e test18.f90 ]
            then
               mv test18.f90 mtest18.f90
            fi ;;
         yes)
            if [ -e mtest18.f90 ]
            then
               mv mtest18.f90 test18.f90
            fi ;;
      esac ;;
esac

case "$PHAML_BLOPEX" in
   no)
      if [ -e test19.f90 ]
      then
         mv test19.f90 mtest19.f90
      fi ;;
   yes)
      if [ -e mtest19.f90 ]
      then
         mv mtest19.f90 test19.f90
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
