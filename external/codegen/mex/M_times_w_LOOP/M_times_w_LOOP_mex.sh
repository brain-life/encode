MATLAB="/N/soft/rhel6/matlab/2015a"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/N/u/ccaiafa/Karst/.matlab/R2015a"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for M_times_w_LOOP" > M_times_w_LOOP_mex.mki
echo "CC=$CC" >> M_times_w_LOOP_mex.mki
echo "CFLAGS=$CFLAGS" >> M_times_w_LOOP_mex.mki
echo "CLIBS=$CLIBS" >> M_times_w_LOOP_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> M_times_w_LOOP_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> M_times_w_LOOP_mex.mki
echo "CXX=$CXX" >> M_times_w_LOOP_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> M_times_w_LOOP_mex.mki
echo "CXXLIBS=$CXXLIBS" >> M_times_w_LOOP_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> M_times_w_LOOP_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> M_times_w_LOOP_mex.mki
echo "LD=$LD" >> M_times_w_LOOP_mex.mki
echo "LDFLAGS=$LDFLAGS" >> M_times_w_LOOP_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> M_times_w_LOOP_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> M_times_w_LOOP_mex.mki
echo "Arch=$Arch" >> M_times_w_LOOP_mex.mki
echo OMPFLAGS= >> M_times_w_LOOP_mex.mki
echo OMPLINKFLAGS= >> M_times_w_LOOP_mex.mki
echo "EMC_COMPILER=gcc" >> M_times_w_LOOP_mex.mki
echo "EMC_CONFIG=optim" >> M_times_w_LOOP_mex.mki
"/N/soft/rhel6/matlab/2015a/bin/glnxa64/gmake" -B -f M_times_w_LOOP_mex.mk
