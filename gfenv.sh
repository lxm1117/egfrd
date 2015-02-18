export PYTHONPATH=$(pwd)
export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH
env | grep PYTHONPATH
env | grep LD_LIBRARY_PATH
