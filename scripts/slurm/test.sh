#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=3:00:00
#SBATCH --mem=32gb
#SBATCH --output=slurm_out/pixel_%A_%a.out
#SBATCH --error=slurm_out/pixel_%A_%a.err

set -x

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

row_index=$1
N=$2
seed=$3

shift 3

param1_min=77.0
param1_max=200.0
param1_step=$(echo "scale=6; ($param1_max - $param1_min)/($N - 1)" | bc -l)

param2_min=0.01
param2_max=0.2
param2_step=$(echo "scale=6; ($param2_max - $param2_min)/($N - 1)" | bc -l)

param1_vals=()
param2_vals=()

for (( i=0; i<${N}; i++ )); do
    v1=$(echo "scale=6; $param1_min + $i * $param1_step" | bc -l)
    v2=$(echo "scale=6; $param2_min + $i * $param2_step" | bc -l)
    param1_vals+=("$v1")
    param2_vals+=("$v2")
done

task_id=${SLURM_ARRAY_TASK_ID}
p1=${param1_vals[$row_index]}
p2=${param2_vals[$task_id]}

file="row_${row_index}_${task_id}"
task_seed="$(( seed + row_index * task_id + N ))"

"$@" \
    --file_name="$file" \
    --seed="$task_seed" \
    --paramName1="temp" \
    --paramName2="sigma" \
    --paramValue1="${p1}" \
    --paramValue2="${p2}" \

    #!/bin/bash
#SBATCH --partition=gpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:A40:1
#SBATCH --time=01:30:00
#SBATCH --mem=16G

module load devel/cuda
#conda activate py311

python3 scripts/python_scripts/surrogate_model.py