#!/bin/bash

execution_successful=false
profile=$1

# If `profile` is empty, set `profile` to "cluster"
if [ -z "$profile" ]; then
    profile="cluster"
fi


# While not execution_successful
while [ $execution_successful == false ]; do
  snakemake --profile $profile
  error_code=$?

  # Execute snakemake again if error code is not 0
  if [ $error_code == 0 ]; then
    execution_successful=true
  fi
done
