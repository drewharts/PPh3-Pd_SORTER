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
	for i in $1/*.opt.mol ; do cat $i ; echo -e "\$\$\$\$" ; done > $2
fi