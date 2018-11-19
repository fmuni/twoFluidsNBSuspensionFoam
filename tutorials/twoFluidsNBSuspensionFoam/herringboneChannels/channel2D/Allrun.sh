#!/bin/sh

currDir=$PWD

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./Allclean.sh

#Run the Case
runApplication blockMesh
#runApplication snappyHexMesh -overwrite
#runApplication transformPoints -scale '(1e-06 1e-06 1e-06)'
runApplication decomposePar
runParallel twoFluidsNBSuspensionFoam
