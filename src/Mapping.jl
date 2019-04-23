module Mapping
using FastaIO

push!(LOAD_PATH,@__DIR__)
using CommonUtils
using JSON

musclepath = fasttreepath = joinpath(@__DIR__,"..","binaries","muscle3.8.31_i86linux64")
if Sys.iswindows()
      musclepath = joinpath(@__DIR__,"..","binaries","muscle3.8.31_i86win32.exe")
end

export unalignedpos
function unalignedpos(sequence, pos)
      len = length(sequence)
      u = 1
      for i=1:len
            if sequence[i] != '-'
                  if i ==  pos
                        return u
                  end
                  u += 1
            end
      end
      return 0
end

export createmapping
function createmapping(fastafile, mappingsequence)
      names = []
      sequences = []
      seqnametoindex  = Dict{AbstractString,Int}()
      FastaIO.FastaReader(fastafile) do fr
            seqindex = 1
            for (desc, seq) in fr
                  len = length(seq)
                  push!(names,desc)
                  push!(sequences, seq)
                  seqnametoindex[desc] = seqindex
                  seqindex += 1
            end
      end

      fastastring = open(fastafile) do file
            read(file, String)
      end
      cachepath = joinpath(@__DIR__,"..","cache")
      mkpath(cachepath)

      key = CommonUtils.sha256base36(string(CommonUtils.sha256base36(fastastring), CommonUtils.sha256base36(mappingsequence)))
      cachefile = joinpath(cachepath, string(key, ".muscle.mapping.json"))
      if isfile(cachefile)
            fin = open(cachefile, "r")
            jsondict = JSON.parse(fin)
            close(fin)
            return convert(Array{Int,1}, jsondict["mapping"]), convert(Array{Int,1}, jsondict["revmapping"])
      else
            alignmentfile1, outfile = mktemp()
            s = 1
            for seq in sequences
                  write(outfile, string(">align1_", s, "\n"))
                  #write(outfile, seq, "\n")
                  write(outfile, replace(seq, "-" => "N"), "\n")
                  s += 1
            end
            close(outfile)

            alignmentfile2, outfile = mktemp()
            write(outfile, ">align2_1\n")
            write(outfile, replace(mappingsequence, "-" => "N"),"\n")
            close(outfile)

            mappingalignmentmuscle, outfile = mktemp()
            close(outfile)
            run(`$(musclepath) -profile -in1 $alignmentfile1 -in2 $alignmentfile2 -out $mappingalignmentmuscle`)
            #rm(alignmentfile1)
            #rm(alignmentfile2)

            align1 = ""
            align2 = ""
            for (desc, seq) in FastaIO.FastaReader(mappingalignmentmuscle)
                  if startswith(desc, "align1_")
                        align1 = seq
                  elseif startswith(desc, "align2_")
                        align2 = seq
                  end
            end
            #rm(mappingalignmentmuscle)

            #mapping = Dict{Int,Int}()
            mapping = zeros(Int, length(sequences[1]))
            revmapping = zeros(Int, length(mappingsequence))
            println(align1)
            println(align2)
            for i=1:length(align1)
                  a1 = unalignedpos(align1,i)
                  a2 = unalignedpos(align2,i)
                  if a1 > 0
                        mapping[a1] = a2
                  end
                  if a2 > 0
                        revmapping[a2] = a1
                  end
            end
            #println(mapping)
            #println(revmapping)
            jsondict = Dict{String,Any}()
            jsondict["mapping"] = mapping
            jsondict["revmapping"] = revmapping
            fout = open(cachefile, "w")
            write(fout, JSON.json(jsondict))
            close(fout)

            return mapping, revmapping
      end
end

export createmapping2
function createmapping2(fastafile, fastafile2)
      names = []
      sequences = []
      seqnametoindex  = Dict{AbstractString,Int}()
      FastaIO.FastaReader(fastafile) do fr
            seqindex = 1
            for (desc, seq) in fr
                  len = length(seq)
                  push!(names,desc)
                  push!(sequences, seq)
                  seqnametoindex[desc] = seqindex
                  seqindex += 1
            end
      end
      alignmentfile1 = string(fastafile,".in1")
      outfile = open(alignmentfile1,"w")
      s = 1
      for seq in sequences
            write(outfile, string(">align1_", s, "\n"))
            write(outfile, replace(seq, "-" => "N"), "\n")
            s += 1
      end
      close(outfile)

      names = []
      sequences = []
      seqnametoindex  = Dict{AbstractString,Int}()
      FastaIO.FastaReader(fastafile2) do fr
            seqindex = 1
            for (desc, seq) in fr
                  len = length(seq)
                  push!(names,desc)
                  push!(sequences, seq)
                  seqnametoindex[desc] = seqindex
                  seqindex += 1
            end
      end
      alignmentfile2 = string(fastafile2,".in2")
      outfile = open(alignmentfile2,"w")
      s = 1
      for seq in sequences
            write(outfile, string(">align2_", s, "\n"))
            write(outfile, replace(seq, "-" => "N"), "\n")
            s += 1
      end
      close(outfile)

      mappingalignmentmuscle= string(fastafile,".musclemapping")
      run(`$(musclepath) -profile -in1 $alignmentfile1 -in2 $alignmentfile2 -out $mappingalignmentmuscle`)

      align1 = ""
      align2 = ""
      for (desc, seq) in FastaIO.FastaReader(mappingalignmentmuscle)
            if startswith(desc, "align1_")
                  align1 = seq
            elseif startswith(desc, "align2_")
                  align2 = seq
            end
      end

      mapping = Dict{Int,Int}()
      revmapping = Dict{Int,Int}()
      for i=1:length(align1)
            a1 = unalignedpos(align1,i)
            a2 = unalignedpos(align2,i)
            if a1 > 0
                  mapping[a1] = a2
            end
            if a2 > 0
                  revmapping[a2] = a1
            end
      end


      println(mapping)
      println(revmapping)

      return mapping, revmapping
end

#=
mapping, revmapping = createmapping2("/media/michael/Sandisk500GB/data/hiv1.fas.norm50.structure.align","/media/michael/Sandisk500GB/data/pNL4-3_trunc.fasta")
for i=1:length(mapping)
println(i,"\t",mapping[i])
end
for i=1:length(revmapping)
println(i,"\t",get(revmapping,i,-1))
end=#

#=
table = readtable("/media/michael/Sandisk500GB/data/NL4-3 SHAPE reactivities.csv")
println(table)
println(table[3])


sequence, paired = readctfile("/media/michael/Sandisk500GB/data/hiv1-SHAPE-MaP.ct")
fastafile = "/media/michael/Sandisk500GB/data/hiv1.fas.norm50.structure.align"
#println(createmapping(fastafile, sequence))=#
end
