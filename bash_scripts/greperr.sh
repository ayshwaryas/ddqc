#!/bin/bash


while getopts j: option
do
case "${option}"
in
j) JOB=${OPTARG};;
esac
done

grep -rn --ignore-case "error\|kill" ~/logs/*.o$JOB* 

