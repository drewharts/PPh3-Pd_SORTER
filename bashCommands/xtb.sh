#!/bin/bash
if [[ $1 == "" ]]
then
 echo "Add the path to the xyz file as the first argument."
else
  for file in $1/*.mol; do xtb $file --opt > temp.out; echo $file | sed 's/.mol/.opt.mol/g'; mv xtbopt.mol `echo $file | sed 's/.mol/.opt.mol/g'`; done
fi
