#!/bin/bash
################ COMPILER OPTIONS (INTEL, ODYSSEY, and EMBRAPA ONLY) #######################
#------------------------------------------------------------------------------------------#
# A/B. Pickiest - Use this whenever you change arguments on functions and subroutines.     #
#                 This will perform the same tests as B but it will also check whether all #
#                 arguments match between subroutine declaration and subroutine calls.     #
#                 WARNING: In order to really check all interfaces you must compile with   #
#                          this option twice:                                              #
#                 1. Compile (./install.sh A)                                              #
#                 2. Prepare second compilation(./2ndcomp.sh)                              #
#                 3. Compile one more time (./install.sh B)                                #
#                 If the compilation fails either at step 3, then your code has interface  #
#                 problems. If it successfully compiles, then the code is fine for         #
#                 interfaces.                                                              #
# C. Pickiest with no interface - This will compile fast but the run will be slow due to   #
#    the -O0 option. However, by setting -O0 you will take full advantage of the intel     #
#    debugger.                                                                             #
#    Ideally, you should compile your code with this option whenever you make any changes. #
#    Note, however, that if you change arguments you should first try A.                   #
# D. Fast check - This will check pretty much the same as C, but it will not set up        #
#    anything for the debugger. Use this only if you really don't want to deal with idb or #
#    if you have a good idea of which problem you are dealing with.                        #
# E. Fast - This is all about performance, use only when you are sure that the model has   #
#           no code problem, and you want results asap. This will not check for any        #
#           problems, which means that this is an option suitable for end users, not de-   #
#           velopers.                                                                      #
#------------------------------------------------------------------------------------------#

#----- Define the number of arguments. ----------------------------------------------------#
nargs=$#
args=$@
#------------------------------------------------------------------------------------------#

# Initialize vars
CLEAN=""
KIND=""
PLATFORM=""
OPT=""
USE_GIT=true

# Argument parsing
while [[ $# > 0 ]]
do
key="$1"
   case $key in
   -p|--platform)
      PLATFORM="$2"
      shift # past argument
      ;;
   -k|--kind)
      KIND="$2"
      shift # past argument
      ;;
   -c|--clean)
      CLEAN="clean"
      ;;
   -g|--gitoff)
      USE_GIT=false
      ;;
   *)
      echo "Unknown key-value argument pair."
      exit 2
      ;;
   esac

   shift # past argument or value
done

# Let the user know if defaults are being used
if [ "${PLATFORM}" == "" ]
then
   echo "No platform specified, defaulting to gfortran."
   PLATFORM="gfortran"
fi

if [ "${KIND}" == "" ]
then  
   echo "No optimization level specified, defaulting to E."
   KIND="E"
fi

#Check that HDF5_HOME is set:
if [ "${HDF5_HOME}" == "" ] && [ "${CLEAN}" != "clean" ]
then
   echo " "
   echo "Your HDF5_HOME variable is not set, attempting to set it using 'find / -name hdf5.mod'."
   HDF5_GUESS=`find / -name hdf5.mod 2>/dev/null`

   if [ "${HDF5_GUESS}" == "" ]
   then
      echo "Failed to set HDF5_HOME. Please set this environment variable. Are you sure HDF5 is installed correctly?"
      exit 1
   else
      # Strip '/bin/h5dump' from path
      # Bash 4.2 and above:
      #HDF5_HOME=${HDF5_GUESS::-11}

      # Below Bash 4.2 compatible:
      HDF5_HOME=`echo ${HDF5_GUESS} | rev | cut -c 9- | rev`
      echo "HDF5_HOME set to '${HDF5_HOME}'"
   fi
fi

# Look for necessary hdf5 libraries:
if [ "${HDF5_LIB_PATH}" == "" ] && [ "${CLEAN}" != "clean" ]
then
   echo " "
   echo "Your HFD5_LIB_PATH variable is not set, attempting to set it using 'find / -name libhdf5_fortran.a'."
   HDF5_LIB_GUESS=`find / -name libhdf5_fortran.a 2>/dev/null`
   
   if [ "${HDF5_LIB_GUESS}" == "" ]
   then
      echo "Failed to set HDF5_LIBS. Please set this environment variable. Are you sure HDF5 is installed correctly?"
      exit 1
   else
      # Below Bash 4.2 compatible:
      HDF5_LIB_PATH=`echo ${HDF5_LIB_GUESS} | rev | cut -c 18- | rev`
      echo "HDF5_LIB_PATH set to '${HDF5_LIB_PATH}'"
   fi
fi

# Set opt and bin
case ${KIND} in
['A','B','C','D']*)
   OPT='dbg'
   ;;
['E']*)
   OPT='opt'
   ;;
*)
   # Default to opt
   echo "Compiler optimization not recognized as opt or dbg."
   exit 1
   ;;
esac

# Tag executables with a git version and branch name if possible.
GIT_EXIST=`git rev-parse --is-inside-work-tree`
if [ ${GIT_EXIST} == "true" -a ${USE_GIT} ]
then
   GIT_TAG=`git branch -v | awk '/\*/ {print "-" $2 "-" $3}'`
   GIT_TAG=`echo ${GIT_TAG} | tr -d '()/[]'`
   echo " "
   echo "Git found, it will be used to tag things."
   echo "To disable revision tagging, use --gitoff or -g."
   echo " "
else
   GIT_TAG=''
fi

BIN=bin-${OPT}-${KIND}${GIT_TAG}

# Move to the binary directory
if [ ! -d "$BIN" ]; then
   mkdir ${BIN}
fi
cd ${BIN}


# Link to makefiles, includes, and shell scripts
ln -sf ../make/*.mk ./
ln -sf ../make/Makefile ./
ln -sf ../make/include.mk.${OPT}.${PLATFORM} ./include.mk
ln -sf ../shell/* ./
touch dependency.mk

#----- Launch the compiler. ---------------------------------------------------------------#
make OPT=${OPT} KIND_COMP=${KIND} ${CLEAN} GIT_TAG=${GIT_TAG} HDF5_HOME=${HDF5_HOME} HDF5_LIB_PATH=${HDF5_LIB_PATH}
make_exit_code=$?
#------------------------------------------------------------------------------------------#

if [ ${make_exit_code} != 0 ]
then
   exit 1
else
   echo "Installation Complete."
   exit 0
fi
