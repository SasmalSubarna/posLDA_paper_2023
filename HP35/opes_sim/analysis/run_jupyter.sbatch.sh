#!/bin/bash
#SBATCH --job-name=jupyter-sbatch
#SBATCH -o logs/jupyter-sbatch-%j
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node 4
#SBATCH --mem=20GB
#SBATCH --time=48:00:00

source ~/.bashrc
module purge
module load python/intel/3.8.6;
#source ~/pyenv/md-py3.8/bin/activate
source /scratch/ss12902/pyenv/bin/activate

port=$(shuf -i 6000-9999 -n 1)

#/usr/bin/ssh -N -f -R $port:localhost:$port log-0
/usr/bin/ssh -N -f -R $port:localhost:$port log-1
/usr/bin/ssh -N -f -R $port:localhost:$port log-2
/usr/bin/ssh -N -f -R $port:localhost:$port log-3


host=$(hostname)

cat<<EOF

Jupyter server is running on: $(hostname)
Job starts at: $(date)

Step 1 :

If you are working in NYU campus, please open an iTerm window, run command

(or log-0)
ssh -L $port:localhost:$port $USER@log-1.hpc.nyu.edu

If you are working off campus, you should already have ssh tunneling setup through HPC bastion host,
so that you can directly login to greene with command

ssh $USER@greene

Please open a terminal window, and run command

ssh -i ~/.ssh/id_nyu -L $port:localhost:$port $USER@greene.hpc.nyu.edu
ssh -i ~/.ssh/id_nyu -L $port:localhost:$port $USER@$host

Step 2:

Keep the iTerm windows in the previouse step open. Now open browser, find the line with

The Jupyter Notebook is running at: $(hostname)

the URL is something: http://localhost:${port}/?token=XXXXXXXX (see your token below)

you should be able to connect to jupyter notebook running remotly on greene compute node with above url

EOF

unset XDG_RUNTIME_DIR
if [ "$SLURM_JOBTMP" != "" ]; then
    export XDG_RUNTIME_DIR=$SLURM_JOBTMP
fi

jupyter notebook --no-browser --port $port --notebook-dir=$(pwd)
