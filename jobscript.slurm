#!/usr/bin/env bash
#
#SBATCH --job-name=ecc
#SBATCH --output=outputfile.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=24:00:00

for p in 0.33 -0.33 -0.67
do
    for centrality in 10-20% 20-30% 30-40% 40-50%
    do
        # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=psin --p=${p} --tau0=1.0 --centrality=${centrality};
        # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=epsn --p=${p} --tau0=1.0 --centrality=${centrality};
        # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=avgevol --p=${p} --tau0=1.0 --centrality=${centrality};
        # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=eccavgevol --p=${p} --tau0=1.0 --centrality=${centrality};
        # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=plotevol --p=${p} --tau0=1.0 --centrality=${centrality};
        # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=jtnpsintau --p=${p} --tau0=1.0 --centrality=${centrality};
        python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=epsnavg --p=${p} --tau0=1.0 --centrality=${centrality};
    done
done

for p in 0.0 -1.0 -2.0 -3.0
do
    for tau0 in 1.0 0.6 0.2
    do
        for centrality in 10-20% 20-30% 30-40% 40-50%
        do
            # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=psin --p=${p} --tau0=${tau0} --centrality=${centrality};
            # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=epsn --p=${p} --tau0=${tau0} --centrality=${centrality};
            # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=avgevol --p=${p} --tau0=${tau0} --centrality=${centrality};
            # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=eccavgevol --p=${p} --tau0=${tau0} --centrality=${centrality};
            # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=plotevol --p=${p} --tau0=${tau0} --centrality=${centrality};
            # python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=jtnpsintau --p=${p} --tau0=${tau0} --centrality=${centrality};
            python3 run.py --NUM_THREADS=${SLURM_CPUS_PER_TASK} --calculation=epsnavg --p=${p} --tau0=${tau0} --centrality=${centrality};
        done
    done
done