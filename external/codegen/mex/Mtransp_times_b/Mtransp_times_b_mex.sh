MATLAB="/Applications/MATLAB_R2014b.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/ccaiafa/.matlab/R2014b"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for Mtransp_times_b" > Mtransp_times_b_mex.mki
echo "CC=$CC" >> Mtransp_times_b_mex.mki
echo "CFLAGS=$CFLAGS" >> Mtransp_times_b_mex.mki
echo "CLIBS=$CLIBS" >> Mtransp_times_b_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> Mtransp_times_b_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> Mtransp_times_b_mex.mki
echo "CXX=$CXX" >> Mtransp_times_b_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> Mtransp_times_b_mex.mki
echo "CXXLIBS=$CXXLIBS" >> Mtransp_times_b_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> Mtransp_times_b_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> Mtransp_times_b_mex.mki
echo "LD=$LD" >> Mtransp_times_b_mex.mki
echo "LDFLAGS=$LDFLAGS" >> Mtransp_times_b_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> Mtransp_times_b_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> Mtransp_times_b_mex.mki
echo "Arch=$Arch" >> Mtransp_times_b_mex.mki
echo OMPFLAGS= >> Mtransp_times_b_mex.mki
echo OMPLINKFLAGS= >> Mtransp_times_b_mex.mki
echo "EMC_COMPILER=Xcode with Clang" >> Mtransp_times_b_mex.mki
echo "EMC_CONFIG=optim" >> Mtransp_times_b_mex.mki
"/Applications/MATLAB_R2014b.app/bin/maci64/gmake" -B -f Mtransp_times_b_mex.mk
