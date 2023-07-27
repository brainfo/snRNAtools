#!/bin/bash -l

#SBATCH -A sens2022004
#SBATCH --partition=core             # long, fast, etc.
##SBATCH --qos=short
#SBATCH --ntasks=1                   # nb of *tasks* to be run in // (usually 1), this task can be multithreaded (see cpus-per-task)
#SBATCH --nodes=1                    # nb of nodes to reserve for each task (usually 1)
#SBATCH --cpus-per-task=8            # nb of cpu (in fact cores) to reserve for each task /!\ job killed if commands below use more cores
#SBATCH --mem=80G                  # amount of RAM to reserve for the tasks /!\ job killed if commands below use more RAM
#SBATCH --time=0-06:00               # maximal wall clock duration (D-HH:MM) /!\ job killed if commands below take more time than reservation
#SBATCH -o /proj/sens2022004/nobackup/wharf/hongki/hongki-sens2022004/human_placenta/pcos/workflow/logs/clean.%A.%a.out   # standard output (STDOUT) redirected to these files (with Job ID and array ID in file names)
#SBATCH -e /proj/sens2022004/nobackup/wharf/hongki/hongki-sens2022004/human_placenta/pcos/workflow/logs/clean.%A.%a.err  # standard error  (STDERR) redirected to these files (with Job ID and array ID in file names)
#/!\ Note that the ./outputs/ dir above needs to exist in the dir where script is submitted **prior** to submitting this script
##SBATCH --array=1-8                # 1-N: clone this script in an array of N tasks: $SLURM_ARRAY_TASK_ID will take the value of 1,2,...,N
#SBATCH --job-name=snakemake_test_qc        # name of the task as displayed in squeue & sacc, also encouraged as srun optional parameter
#SBATCH --mail-type END              # when to send an email notiification (END = when the whole sbatch array is finished)
#SBATCH --mail-user scilavisher@gmail.com
# ml bioinfo-tools
module purge
ml conda
export CONDA_ENVS_PATH=/proj/sens2022004/nobackup/wharf/hongki/hongki-sens2022004/conda_envs/
# ml snakemake
conda activate scanpy
module purge
snake_base_dir="/proj/sens2022004/nobackup/wharf/hongki/hongki-sens2022004/human_placenta/pcos/workflow"
cd ${snake_base_dir}
# snakemake -s Snakefile -j 1 --cluster-config cluster.yaml --cluster "sbatch SBATCH -A sens2022004 -t {cluster.time} -p {cluster.partition} -n {cluster.n} --mem {cluster.mem} --cpus-per-task {cluster.cpu}"
snakemake -s Snakefile -j 1