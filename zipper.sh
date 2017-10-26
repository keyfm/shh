#!/bin/bash

if [ $# -le 1 ]
then
	echo "Usage: script zips/unzips files in provided folder using parallel"
	exit 1
fi

# Read input from CL
echo $#
switch="off"
while [ $# != "0" ] ; do
 if [[ ${switch} == 'off' ]]; then
  input=""
  zipcmd=${1} # <zip> <unzip>
  switch="on"
 else
  input="$input $1"
 fi
 shift
done

echo "Hallo"
echo $zipcmd
echo $input
echo "BYE"

# Unzip
if [[ ${zipcmd} == 'unzip' ]]; then 
    parallel --jobs 30 "echo -e '\n'{} processing 1>&2 ; bgzip -d {} " ::: ${input}
fi

# zip
if [[ ${zipcmd} == 'zip' ]]; then 
    parallel --jobs 30 "echo -e '\n'{} processing 1>&2 ; bgzip {} " ::: ${input}
fi
    
    

