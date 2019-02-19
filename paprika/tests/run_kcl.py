import os as os
import glob as glob

directories = glob.glob("tmp/windows/*")

for directory in sorted(directories):
    print(directory)
    if os.path.exists(f"{directory}/production.nc"):
        print("...found simulation and skipping")
        continue

    tscc = f"""
#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=2 -q home-gibbs
#PBS -j oe -r n
#PBS -N {directory.split('/')[-1]}

source /home/davids4/amber18_gnu.sh

SCRDIR=/oasis/tscc/scratch/davids4/k-cl/{directory.split('/')[-3]}-{directory.split('/')[-1]}
mkdir -p $SCRDIR

# Need the `-L` to resolve any links.
rsync -avL $PBS_O_WORKDIR/ $SCRDIR/

cd $SCRDIR
pmemd.cuda -O -p k-cl-sol.prmtop -ref k-cl-sol.rst7 -c minimize.rst7 -i production.in -o production.out -r production.rst7 -x production.nc -inf production.mdinfo -e production.mden

"""
    with open(os.path.join(directory, "tscc-driver.sh"), "w") as f:
        f.write(tscc)

    print(f"Running qsub tscc-driver.sh in {directory}")
    sp.call("qsub tscc-driver.sh", cwd=directory, shell=True)

