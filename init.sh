#!/usr/bin/env bash

# Initialization script
## Will only run from git directory 
## MUST BE RERUN IF MOVING GIT DIRECTORY

if [ ! -f "init.sh" ]; then
  echo "ERROR : Must be run from where git directory was cloned"
  exit -1
fi

current_dir=$(pwd)
dock_dir=${current_dir}/bin


for f in ${dock_dir}/*.csh; do
  sed -i "s|REPLACE_ME|$dock_dir|" ${f}
done
