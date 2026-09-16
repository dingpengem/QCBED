#!/bin/bash

#SBATCH --partition=CPUQ
#SBATCH --time=0-4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --job-name="padqcbed"
#SBATCH --output=padqcbed
#SBATCH --mem=4000
#SBATCH --mail-user=elisathr@ntnu.no
#SBATCH --account=share-nv-fys-tem 


module load intel

touch REFINEMENTSUMMARY.DAT
for i in `seq 1 5`
do
	./qcbedms
	cat PARS.DAT REFINEMENTSUMMARY.DAT >> TEMP.DAT
	mv TEMP.DAT REFINEMENTSUMMARY.DAT
	cp out.dat bestout.dat
	./distcorr
	cp parameters.dat bestparameters.dat
	cp in.dat bestin.dat
done
