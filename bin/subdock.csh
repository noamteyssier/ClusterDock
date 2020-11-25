#!/bin/csh -f

if ( $#argv != 1 ) then
    echo
    echo "Submit jobs to SGE as array jobs using specified DOCK version."
    echo
    echo "usage: subdock.csh path/to/dock_executable"
    echo
    exit 1
endif

set dock = "$1"

if ( ! -e dirlist ) then
    echo "Error: Cannot find dirlist, the list of subdirectories!"
    echo "Exiting!"
    exit 1
endif

set dirnum=`cat dirlist | wc -l`
set dock_dir = "/wynton/home/rotation/nteyssier/ClusterDock/bin"

qsub -t 1-$dirnum  ${dock_dir}/rundock.csh "$dock"
