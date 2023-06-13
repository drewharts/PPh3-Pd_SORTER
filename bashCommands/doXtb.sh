#!/bin/bash
if [[ $1 == "" ]]
then
	echo "Add the path to the xyz file as the first argument."
else
  for file in $1/*.xyz; do xtb $file --opt; mv xtbopt.xyz `echo $file | sed 's/.xyz/.opt.xyz/g'`; done
fi