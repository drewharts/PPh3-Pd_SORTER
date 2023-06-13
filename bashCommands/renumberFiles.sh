#!/bin/bash
var=0
if [[ $1 == "" ]]
then
	echo "Add the path to the xyz file as the first argument."
else
  for file in $1/*.xyz; do
    var=$((var+1)); mv "$file" "$1/file_$var.xyz";
    done;
fi
