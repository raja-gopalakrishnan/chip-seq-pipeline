#!/bin/bash

#SBATCH -p priority
#SBATCH -t 18:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -n 1
#SBATCH -e superjob.err
#SBATCH -o superjob.log
#SBATCH -J ChIPSeq-snakemake

#snakemake -n -R `snakemake --lc --li --lp` --latency-wait 300 --rerun-incomplete --cluster-config cluster.yaml --use-conda --jobs 999 --cluster "sbatch -p {cluster.queue} -n {cluster.n} -t {cluster.time} --mem-per-cpu={cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log}"

snakemake -p --latency-wait 300 --rerun-incomplete --cluster-config cluster.yaml --use-conda --jobs 999 --cluster "sbatch -p {cluster.queue} -n {cluster.n} -t {cluster.time} --mem-per-cpu={cluster.mem} -J {cluster.name} -e {cluster.err} -o {cluster.log}"

#"sbatch -p short -n 4 -t 6:00:00 --mem-per-cpu=4G -e snakemake.err -o snakemake.log"
