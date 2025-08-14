#!/bin/bash

# Define parameters
jobscript=~/jupyterlab.slurm
hostAddress=helix.bwservices.uni-heidelberg.de

# Run job
job_id=$(sbatch $jobscript | awk '{print $4}')
echo "jobid: $job_id"

# Outfile name
slurm_out=slurm-${job_id}.out

# Wait for output file
while [ ! -f $slurm_out ]; do   
    sleep 2; 
done

# Wait until url is written in output file
while [ -z ${url} ]; do   
    sleep 1; 
    url=$(grep -o 'http[^ ]*' $slurm_out | head -n 1); 
    done

# Extract hostID and port from output. The pattern assumes a node name with a length of 6 characters and a port with a length of 3, 4 or 5 numbers.
url_pattern="http://([a-z0-9]{6}):([0-9]{3,5})/lab"
if [[ $url =~ $url_pattern ]]; then 
    hostID=${BASH_REMATCH[1]}
    port=${BASH_REMATCH[2]}
    echo "To connect with the JupyterLab kernel, please enter the following into your local commandline: "
    echo "ssh -N -L $port:$hostID:$port ${USER}@$hostAddress"; 
    echo ""
    echo "Note: It is normal that the ssh command doesn't end after providing the credentials. Ending the command would mean ending the local connection to the kernel."
    echo ""
    echo "Afterwards, you can use the URL"
    echo "  http://127.0.0.1:${port}/lab "
    echo ""
    echo "to:"
    echo "- use the kernel in VSCode ('Existing Jupyter Server...', enter URL, enter password, confirm '127.0.0.1', choose kernel) or "
    echo "- open JupyterLab in your browser with the URL"
else
    echo "The needed information couldn't be found in the slurm output. Please contact your support unit if you need help with fixing this problem."
fi
# rm $slurm_out