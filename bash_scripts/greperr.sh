#!/bin/bash

cd ~/method_comparison

while getopts j: option
do
case "${option}"
in
j) JOB=${OPTARG};;
esac
done

grep -rn "Error*" logs/*.o$JOB* 

