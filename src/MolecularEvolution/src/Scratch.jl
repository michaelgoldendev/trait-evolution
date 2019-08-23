push!(LOAD_PATH,string(@__DIR__))

using MolecularEvolution
using FastaIO

function gettreenodes(infile)
    outfile = string(infile,".tsv")
    fout = open(outfile,"w")
    newick = read(open(infile,"r"), String)
    root = gettreefromnewick(newick)
    nodelist = getnodelist(root)
end


namelabeldict = Dict{String,Int}()
nameseqdict = Dict{String,String}()

files = AbstractString[]
#push!(files, "/media/michael/Sandisk500GB/Dropbox/dev/SFTSV/Death_Survive_Seqs_L_P.fasta")
#push!(files, "/media/michael/Sandisk500GB/Dropbox/dev/SFTSV/Death_Survive_Seqs_M_P.fasta")
#push!(files, "/media/michael/Sandisk500GB/Dropbox/dev/SFTSV/Death_Survive_Seqs_S_P.fasta")
push!(files, "D:/Dropbox/dev/SFTSV/Death_Survive_Seqs_L_P.fasta")
#push!(files, "D:/Dropbox/dev/SFTSV/Death_Survive_Seqs_M_P.fasta")
#push!(files, "D:/Dropbox/dev/SFTSV/Death_Survive_Seqs_S_P.fasta")

for fastafile in files
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            tag = match(r"(.+_.+)_.+_.+_.+", desc)[1][2:end]
            if occursin("_D_", desc)
                namelabeldict[tag] = 1
            elseif occursin("_S_", desc)
                namelabeldict[tag] = 0
            end
            nameseqdict[tag] = seq
            push!(sequences, seq)
        end
    end
end

#infile1 = "/media/michael/Sandisk500GB/Dropbox/dev/SFTSV/sftsvunder60s/SFTSV_L_sub60s.muscle.fas.midpoint.nwk"
infile1 = "D:/Dropbox/dev/SFTSV/sftsvunder60s/SFTSV_L_sub60s.muscle.fas.midpoint.nwk"
fout = open(string(infile1,".tsv"),"w")
fout2 = open(string(infile1,".fas"),"w")
nodelist = gettreenodes(infile1)
for node in nodelist
    if isleafnode(node)
        tag = replace(split(node.name, "|")[1], "-" => "_")[2:end]
        #println(tag,"\t",get(namelabeldict, tag, -1))
        if get(namelabeldict, tag, -1) >= 0
            println(fout,node.name, "\t", get(namelabeldict, tag, 0))
            if get(namelabeldict, tag, -1) == 0
                println(fout2, string(">", node.name,"_S_"))
            else
                println(fout2, string(">", node.name,"_D_"))
            end
            println(fout2, nameseqdict[tag])
        end
        #println(tag)
        #=
        if occursin("_D_", node.name)
            println(fout, node.name,"\t",1)
        elseif occursin("_S_", node.name)
            println(fout, node.name,"\t",0)
        end=#
    end
end
close(fout)
close(fout2)
#=
infile1 = "/media/michael/Sandisk500GB/Dropbox/dev/SFTSV/RAxML_SFTSV_L.tree.midpoint.nwk"
nodelist = gettreenodes(infile1)
for node in nodelist
    if isleafnode(node)
        if occursin("_D_", node.name)
            println(fout, node.name,"\t",1)
        elseif occursin("_S_", node.name)
            println(fout, node.name,"\t",0)
        end
    end
end
close(fout)
=#
