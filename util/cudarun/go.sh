#-------------------------------------------------------
# Script to rsync code changes to physlogin, compile
# vampire, copy inputs to compute node, run the job and
# print results to screen
#-------------------------------------------------------

# setup
USER=rfle500
NODE=physnode33
SERV=physlogin
VDIR=vampire-cuda

# tell bash to treat hidden files and directories as normal
shopt -s dotglob

# syncrhonise current directory with target directory
echo "-----------------------------------------------"
echo " Synchronising files to $SERV:$VDIR"
echo "-----------------------------------------------"
RSYNC_RSH="ssh -q" rsync -rcih * $USER@$SERV:$VDIR/

# compile the code on the headnode
echo "-----------------------------------------------"
echo " Compiling on $SERV"
echo "-----------------------------------------------"
ssh -qt $USER@$SERV $CMD /bin/bash -l $VDIR/util/cudarun/vcompile.sh -d $VDIR

echo "-----------------------------------------------"
echo " Executing on $NODE"
echo "-----------------------------------------------"
# run the code on the compute node -l login shell -q quiet mode to supress login message -t for terminal
ssh -qt $USER@$SERV /bin/bash -l $VDIR/util/cudarun/vrun.sh -u $USER -d $VDIR -n $NODE
