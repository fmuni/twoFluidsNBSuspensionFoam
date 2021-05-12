#!/bin/sh

currDir=$PWD

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./Allclean.sh

cp -r 0.init 0

#Run the Case
runApplication blockMesh
runApplication setFields
runApplication twoFluidsNBSuspensionFoam
