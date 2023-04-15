#!/usr/bin/bash
#SBATCH --output=log_%x_%j
#SBATCH --partition=YOUR_PARTITION
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
## these simulations, especially the postprocessing,
## take a lot of memory per processor, so you might
## have to use less processors than those available
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --time=20:59:00
#SBATCH --mail-type=ALL

# here we load libraries and python packages
# these are installation dependent
# gsl
module load use.own gsl2/2.7
# pythia 8.303
source /home/hireaction/inghirami/enable_pythia8303
# Anaconda python distribution, version 2021.05
source /home/hireaction/inghirami/enable_conda

export LC_NUMERIC="en_US.UTF-8"

idrun=Au_Elab_1_23_centr_$SLURM_JOB_ID
# path to the SMASH compiled binary
sd=/scratch/hireaction/inghirami/hydro_conditions/smash_bin
# name of the path executable
smash_exe=$sd/smash_fuchs
# config file
smash_conf=$sd/config_Au_Elab_1_23_centr_ncl.yaml
# name of the python script to postproces a bunch of events
analysis_script=analyze_data_ncl.py
# name of the python script to combine the results of the postprocessing
# of various bunches of events
combine_script=combine_data.py

cd $sd

mkdir -p $idrun\_tmp_results

for kk in $(seq 0 9)
do
  echo Iteration $kk

  for i in $(seq 1 $SLURM_NTASKS)
  do
  outdir=out$idrun\_$i
  $smash_exe -i $smash_conf -o $outdir > log_$outdir&
  done
  wait
  sleep 25

  for i in $(seq 1 $SLURM_NTASKS)
  do
  outdir=out$idrun\_$i
  python3 $analysis_script $outdir $idrun\_tmp_results/run$i &
  done
  wait

  python3 $combine_script $idrun\_tmp_results/run* $idrun\_tmp_results/results_step$kk
  sleep 25

  for i in $(seq 1 $SLURM_NTASKS)
  do
  outdir=out$idrun\_$i
  rm -rf $outdir $idrun\_tmp_results/run$i
  done
done
python3 $combine_script $idrun\_tmp_results/results_step* results_$idrun
wait


