# Trait association model

### Installing and running the model:
1. [Download](https://julialang.org/downloads/) and install Julia 1.0 or higher.
2. [Download](https://github.com/michaelgoldendev/trait-evolution/archive/master.zip)'s the source code and extract it to a folder on your computer.
3. Navigate to the source folder using a Windows command prompt or a Linux/Mac terminal:
```
cd /path_to_code/trait-evolution/src/
```
4. Run the model using the path to your julia executable, the path of your alignment file, the path of your tree file, and a comma-seperated list including the primary and secondary annotations:
```
/path_to_julia/bin/julia ParallelEvolution.jl --mlmodel --alignment ../data/H7NX/HA/H7_HA_alnP.fas --tree ../data/H7NX/HA/H7_HA_alnP.fas.nwk --annotations HP,LP
```

The order of the arguments in the annotations list is important, the primary annotation should be listed first, followed by the secondary annotation. Each sequence in the specified alignment and corresponding tree should be tagged with an annotation by simply including the annotation in the sequence name. In the alignment used above (`../data/H7NX/HA/H7_HA_alnP.fas`) each sequence is tagged with either HP (the primary annotation) or LP (the secondary annotation).

By default the model will write the results to a .csv file in the same folder as your alignment (`../data/H7NX/HA/` in the example above).
