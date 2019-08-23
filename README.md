# Trait association model

### Installing and running the model:
1. [Download](https://julialang.org/downloads/) and install Julia 1.0 or higher, if it isn't already installed on your system.
2. [Download](https://github.com/michaelgoldendev/trait-evolution/archive/master.zip)'s this project's source code and extract it to a folder on your computer.
3. Navigate to the source folder using a Windows command prompt or a Linux/Mac terminal:
```
cd /path_to_messi/trait-evolution/src/
```
4. Run the model using the path to your julia executable and the path of your alignment file:
```
/path_to_julia/bin/julia ParallelEvolution.jl --mlmodel --alignment ../data/H7NX/HA/H7_HA_alnP.fas --tree ../data/H7NX/HA/H7_HA_alnP.fas.nwk --annotations LP,HP
```
(the above example should work out of the box)

By default the model will write your results to a .csv file in the same folder as your alignment (../data/H7NX/HA/ in the example above).
