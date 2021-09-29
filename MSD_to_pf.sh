#!/bin/bash

#Dependencias

echo "Searching dependencies ..."
if command -v gnuplot >/dev/null 2>&1 ; then
    echo "gnuplot...OK"
else
    echo "gnuplot...NO"
    echo "Installing gnuplot"
    sudo apt-get -y install gnuplot

fi

# this script is for calculating intrinsic permeability (pf) from MSD(t) obtained from watRpore.R.
#Calculations are made for a tetrameric assembly, so a pf will be computed for each monomer in each trajectory.

echo "Running........"
for j in A B C D
do
echo -e "f(x)=a*x+b
fit f(x) './MSD_`echo $j`.txt' via a,b" > slope.txt

echo -e "\n############# $(date) #############" >> Data.txt
echo -e "######################## CHAIN ${j} #######################\n" >> Data.txt
gnuplot slope.txt &>> Data_raw
awk '{print $3*2.99003322259e-23/2e-12}' fit.log | grep -vE "0$">> Data.txt

# VERY IMPORTANT: the timestep of each frame should be taken in account when calculating pf.
#In this example trajectory was written every 1 ps, so the molar volume of water is divided by 2e-12.
#If every frame would represent 10 ps, then you will need to divide by 2e-11. 
#In this way pf units will be cm^3*s^-1 Refer to Zhu et al. 2004.

rm -rf fit.log slope.txt
done

if [ $? -eq 0 ]
	then
		echo -e "Script done without error \U1F603\n"
	else
		echo -e  "Error \U1F615\n"
fi


echo -e " \n#########################  RAW  ########################\n" >> Data.txt
cat Data_raw >> Data.txt && rm Data_raw
