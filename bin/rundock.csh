#$ -S /bin/csh
#$ -cwd
#$ -j yes
#$ -o rundock.log
#$ -v DOCKBASE
#$ -q long.q


if ( $#argv != 1 ) then
    echo
    echo "Run a dock job inside SGE."
    echo
    echo "usage: rundock.csh path/to/dock_executable"
    echo
    exit 1
endif

set dock = "$1"
set dock_dir = "/wynton/home/rotation/nteyssier/ClusterDock/bin"
#set localdock = `basename $dock`

#echo "Copying DOCK to directory: $PWD"
${dock_dir}/get_dock_files.csh

# switch to subdirectory
echo "Starting dock in directory: $PWD"
set dirarray=`cat dirlist`
if ! $?SGE_TASK_ID then
    set SGE_TASK_ID=1
endif
set pth=$dirarray[$SGE_TASK_ID]
cd $pth

echo "HOST: "`hostname` > stderr
# get the real path with this pwd madness (i.e. dock68 instead of dockenv) 
pushd "$dock:h" > /dev/null
set real_dir = `pwd`
popd > /dev/null
echo "DOCK: $real_dir/$dock:t" >> stderr
if ( $?JOB_ID ) then 
    echo "JOB ID: $JOB_ID" >> stderr
endif
if ( $?SGE_TASK_ID) then
    echo "TASK ID: $SGE_TASK_ID" >> stderr
endif
# now run dock

$dock INDOCK >>& stderr
exit $status
