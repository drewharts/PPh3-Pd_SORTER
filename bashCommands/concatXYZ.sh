#!/bin/bash
if [[ $1 == "" ]]
then
	echo "Add the path to the xyz files as the first argument."
	echo "Add the name of the concatenated file as the second argument."
elif [[ $2 == "" ]]
then
	echo "Add the path to the xyz files as the first argument."
	echo "Add the name of the concatenated file as the second argument."
else
	for i in $1/*.xyz ; do cat $i ; echo -e "\n" ; done > $2
fi
