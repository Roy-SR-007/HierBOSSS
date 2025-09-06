#!/bin/bash
#SBATCH --time=0-01:00:00           # Request 1 hour of wall time
#SBATCH --nodes=1		    # Request 1 node
#SBATCH --ntasks=1                  # Only open 1 instance of the server
#SBATCH --cpus-per-task=4           # Use 4 cores
#SBATCH --mem=8G                    # Use 8G of RAM
#SBATCH --partition=short	    # Use short partition (adjust if needed for longer)
##SBATCH --gres=gpu:1               # Use 1 GPU - the double hash means this is not active 
#SBATCH --output=server.output.%j   # Send output to file named 'server.output.{jobid}'

module purge
module load RStudio-Server

# Use a pre-built script to launch RStudio
/modules/containers/rstudio/start_rstudio_4.3.2.sh

# The log-in credentials can be found in the file 'server.output.{jobid}'
# after the job launches, e.g. 'cat server.output.232', where 232 is your job id.