#!/bin/csh -f

if ( $#argv != 0 ) then
    echo
    echo "Submit jobs to SGE as array jobs."
    echo
    echo "usage: submit.csh"
    echo
    exit 1
endif

set dock_dir = "/wynton/home/rotation/nteyssier/ClusterDock/bin"

${dock_dir}/subdock.csh ${dock_dir}/dock.csh
exit $status
