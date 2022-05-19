#!/bin/bash

cd complex

for w in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
  cd $w

    sbatch <<-_EOF
#!/bin/bash
#SBATCH --job-name=TIcpx${w}
#SBATCH --output=output${w}.txt
#SBATCH --ntasks=8
#SBATCH --time=10-00:00:00
#SBATCH -p bolo,labo1,toul,lisboa,achil
#SBATCH -N 1


module purge
module load amber/16

mpirun -np 8 pmemd.MPI -i min.in -c ti.rst7 -ref ti.rst7 -p ti.parm7 \
  -O -o min.out -inf min.info -e min.en -r min.rst7 -l min.log

mpirun -np 8 pmemd.MPI -i heat.in -c min.rst7 -ref ti.rst7 -p ti.parm7 \
  -O -o heat.out -inf heat.info -e heat.en -r heat.rst7 -x heat.nc -l heat.log

mpirun -np 8 pmemd.MPI -i ti.in -c heat.rst7 -p ti.parm7 \
  -O -o ti001.out -inf ti001.info -e ti001.en -r ti001.rst7 -x ti001.nc \
     -l ti001.log

_EOF

  cd ..
done
