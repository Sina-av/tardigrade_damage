# script to update the file from the files inside the sandbox

echo "UPDATING MATERIAL FILES ..." 
scp ../tardigrade_damage.sbx/projects/tardigrade/include/materials/G*.h ./materials/
scp ../tardigrade_damage.sbx/projects/tardigrade/src/materials/G*.C ./materials/
echo "DONE"
echo ""
echo "UPDATING KERNEL FILES ..."
scp ../tardigrade_damage.sbx/projects/tardigrade/include/kernels/G*.h ./kernels/
scp ../tardigrade_damage.sbx/projects/tardigrade/src/kernels/G*.C ./kernels/
echo "DONE"
