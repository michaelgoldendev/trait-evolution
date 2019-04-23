using FastaIO

fastafile = "H7_alnP.fasta"
sequences = AbstractString[]
names = AbstractString[]
FastaIO.FastaReader(fastafile) do fr
    for (desc, seq) in fr
        push!(names,desc)
        push!(sequences, seq)
    end
end

infile = open("H7_RAxML_CS.taxlabels", "r")
namecolordict = Dict{AbstractString,Int}()
for line in readlines(infile)
    m = match(r"'(.+)'(\[.+(#.+)\])", strip(line))
    #println(m[1],"\t",m[3])
    if m[3] == "#000000"
        namecolordict[m[1]] = 1
    else
        namecolordict[m[1]] = 2
    end
end
close(infile)
outfile = open(string(fastafile,".norm.fas"), "w")
for (name,sequence) in zip(names,sequences)
    newname = replace(name, "/", ".")
    if namecolordict[name] == 1
        newname = string(newname, ".LP.x")
    elseif  namecolordict[name] == 2
        newname = string(newname, ".HP.x")
    end
    println(outfile, ">", newname)
    println(outfile, sequence)
end
close(outfile)
