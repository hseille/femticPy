#! /bin/bash

project=$1
cd $project
cp input_data/{topography.dat,bathymetry.dat,coast_line.dat} meshGen/files/
cd meshGen


#-----------------------"
#-----   STEP 1    -----"
#-----------------------"

echo "STEP 1... "
cp files/{control.dat,coast_line.dat,analysis_domain.dat,observing_site.dat} .
/scratch3/sei029/femtic/makeTetraMesh -stp 1
echo "   done!"
sleep 1

#-----------------------"
#-----   STEP 2    -----"
#-----------------------"

echo "STEP 2... "
/scratch3/sei029/femtic/makeTetraMesh -stp 2
echo "   done!"
sleep 1

#-----------------------"
#-----   STEP 3    -----"
#-----------------------"

echo "STEP 3... "
cp files/{topography.dat,bathymetry.dat} .
/scratch3/sei029/femtic/makeTetraMesh -stp 3
echo "   done!"
sleep 1

#-----------------------"
#-----   STEP 4    -----"
#-----------------------"

echo "STEP 4... "
/scratch3/sei029/femtic/makeTetraMesh -stp 4
sed -i "$(( $(wc -l < output.poly) ))s/.*/2/" output.poly
echo "1   0.0 0.0 -40.0  10 1e9" >> output.poly
echo "2  0.0 0.0  10.0  20 1e9" >> output.poly
#echo "3   0.0 0.0  40.0  30 1e9" >> output.poly
echo "   done!"
sleep 1


#-----------------------"
#-----   STEP 5    -----"
#-----------------------"

echo "STEP 5... "
/scratch3/sei029/femtic/tetgen -nVpYAakq3.0/0 output.poly
echo "   done!"
sleep 1


echo "fix output.1.ele if necessary... "
module load python
python fix_output1.py
sleep 1
echo "   done!"



#-----------------------"
#-----   STEP 6    -----"
#-----------------------"

echo "STEP 6... "
cp files/{obs_site.dat,makeMtr.param} .
for i in 1 2 3 4 5
do
/scratch3/sei029/femtic/makeMtr output.$i
/scratch3/sei029/femtic/tetgen -nmpYVrAakq3.0/0 output.$i
done
echo "   done!"
sleep 1


#-----------------------"
#-----   STEP 7    -----"
#-----------------------"

echo "STEP 7... "
cp files/resistivity_attr.dat .
/scratch3/sei029/femtic/TetGen2Femtic  output.6
echo "   done!"


cp {resistivity_block_iter0.dat,mesh.dat} ../inversion/

