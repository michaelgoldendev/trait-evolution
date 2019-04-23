push!(LOAD_PATH,"../juliamolev/src/")
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using CTMCs
using LG

treefile = "H7_Genome/M/H7_M1_alnP.fasta.nwk"
newickin = open(treefile,"r")
newickstring = strip(readlines(newickin)[1])
close(newickin)
doroottree = true
root = gettreefromnewick(newickstring)
if doroottree
     root = roottree(gettreefromnewick(newickstring), 1)
end
nodelist = getnodelist(root)
for node in nodelist
    if node.branchlength <= 1e-4
        node.branchlength = 1e-4
    end
end
seqnametonodeindex = Dict{AbstractString,Int}()
nodeindex = 1
for node in nodelist
    node.nodeindex = nodeindex
    if node.name != ""
        seqnametonodeindex[node.name] = node.nodeindex
    else
        node.name = string("[",node.nodeindex,"]")
    end
    nodeindex += 1
end
#=
seqindextonodeindex = zeros(Int,length(sequences))
seqindex = 1
for (taxon,sequence) in zip(names, sequences)
    seqindextonodeindex[seqindex] = seqnametonodeindex[taxon]
    nodelist[seqnametonodeindex[taxon]].seqindex = seqindex
    seqindex += 1
end=#


rng = MersenneTwister(10498012421321)
paths,times = modifiedrejectionsampling(rng, LGmatrix, 1, 20, nothing)
println(paths,times)
