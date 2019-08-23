mutable struct Alignment
    sequences::Array{AbstractString,1}
    taxa::Array{AbstractString,1}
    charactertype::Int
    noncoding::Array{Int,1}
    coding1::Array{Int,1}
    coding2::Array{Int,1}
    coding3::Array{Int,1}
    codonindices::Array{Int,1}
    treestrings::Array{AbstractString,1}

    function Alignment()
        new(AbstractString[],AbstractString[],0, Int[], Int[], Int[], Int[], Int[],AbstractString[])
    end
end

export getalignmentfromfastaandnewick
function getalignmentfromfastaandnewick(fastafile::AbstractString, newickfile::AbstractString, startpos::Int=0, endpos::Int=0)
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            if startpos <= 0 || endpos <= 0
                push!(sequences, seq)
            else
                push!(sequences, seq[startpos:endpos])
            end
        end
    end

    newickin = open(newickfile,"r")
    newickstring = strip(readlines(newickin)[1])
    close(newickin)

    alignment = Alignment()
    alignment.sequences = sequences
    alignment.taxa = names
    alignment.treestrings = AbstractString[newickstring]
    datalen = length(alignment.sequences[1])
    if length(alignment.noncoding) == 0 && length(alignment.coding1) == 0 && length(alignment.coding2) == 0 && length(alignment.coding3) == 0
        for i=1:datalen
            push!(alignment.noncoding, i)
        end
    end
    alignment.codonindices = getcodonindices(datalen, alignment.noncoding, alignment.coding1, alignment.coding2, alignment.coding3)
    return alignment
end

function savefasta(alignment::Alignment, outfile::AbstractString)
    out = open(outfile, "w")
    for (taxon, seq) in zip(alignment.taxa, alignment.sequences)
        write(out, string(">", taxon,"\n"))
        write(out, string(seq,"\n"))
    end
    close(out)
end

function savenexus(alignment::Alignment, outfile::AbstractString)
    out = open(outfile, "w")
    write(out, "#NEXUS\n\n")

    write(out, "BEGIN TAXA;\n")
    write(out, string("\tDIMENSIONS NTAX = ", string(length(alignment.taxa)), ";\n"))
    write(out, "\tTAXLABELS\n")
    write(out, string("\t\t", join([string("'",t,"'") for t in alignment.taxa], " "), ";\n"))
    write(out, "END;\n\n")

    maxtaxalength = 0
    for taxon in alignment.taxa
        maxtaxalength = max(maxtaxalength, length(taxon))
    end
    maxtaxalength += 5

    write(out, "BEGIN MATRIX;\n")
    for (name, seq) in zip(alignment.taxa, alignment.sequences)
        write(out, string("\t", rpad(string("'",name,"'"),maxtaxalength) , "\t",seq, "\n"))
    end
    write(out, ";\n")
    write(out, "END;\n\n")

    write(out, "BEGIN TREES;\n")
    treeindex = 1
    for treestring in alignment.treestrings
        write(out, string("\tTREE tree", treeindex," = ", treestring,";\n"))
        treeindex += 1
    end
    write(out, "END;\n\n")

    write(out, "BEGIN CODONS;\n")
    write(out, "\tCODONPOSSET codonpositions = \n")
    write(out, string("\tN: ", join(alignment.noncoding, ","),"\n"))
    write(out, string("\t1: ", join(alignment.coding1, ","),"\n"))
    write(out, string("\t2: ", join(alignment.coding2, ","),"\n"))
    write(out, string("\t3: ", join(alignment.coding3, ","),"\n"))
    write(out, "\t;\n")
    write(out, "END;\n\n")

    close(out)
end

function loadnexus(infile::AbstractString)
    fin = open(infile, "r")
    blocks = AbstractString[]

    alignment = Alignment()
    spl = split(readstring(fin), r"(\s)*END;")
    for s in spl
        regexbegin = match(r".*BEGIN(\s)+([^;]*);([\s\S]*)", normalize_string(s))
        if regexbegin != nothing
            blockname = regexbegin[2]
            content = regexbegin[3]
            #=
            println("-------------------------------------------------------------------")
            println("B:",blockname)
            println("-------------------------------------------------------------------")
            println(content)
            println("-------------------------------------------------------------------")
            =#
            if uppercase(blockname) == "MATRIX" ||  uppercase(blockname) == "CHARACTERS"
                seqdict = Dict{AbstractString,AbstractString}()
                splcontents = split(strip(content),r"[\r\n]+")
                startseqs = uppercase(blockname) == "MATRIX"
                for splcontent in splcontents
                    #println("%%%%",splcontent)
                    splsplcontent = split(splcontent)
                    if startseqs && length(splsplcontent) >= 2
                        taxon = strip(splsplcontent[1], '\'')
                        seqdict[taxon] = string(get(seqdict, taxon, ""), strip(splsplcontent[2],[';']))
                        push!(alignment.taxa, taxon)
                    end
                    if startswith(splcontent,"MATRIX")
                        startseqs = true
                    end
                end
                for taxa in alignment.taxa
                    push!(alignment.sequences, seqdict[taxa])
                end
            elseif uppercase(blockname) == "TREES"
                splcontents = split(strip(content),r"[\r\n]+")
                for splcontent in splcontents
                    if startswith(uppercase(strip(splcontent)), "TREE")
                        regextree = match(r".*=\s*(.*;?).*;", splcontent)
                        push!(alignment.treestrings, regextree[1])
                    end
                end
            elseif uppercase(blockname) == "CODONS"
                splcontents = split(strip(content),"\n")
                for splcontent in splcontents
                    posline = strip(splcontent)
                    regexcodon = match(r"(.:)\s*(.*)\s*;*", posline)
                    if regexcodon != nothing
                        postype = regexcodon[1]
                        posspl = split(regexcodon[2], ",")
                        for posstr in posspl
                            if contains(posstr, "-")
                                regexposstr = match(r"([0-9]+)-([0-9]+)(\\3)?", posstr)
                                modulo3 = regexposstr[3] != nothing
                                startpos = parse(Int,regexposstr[1])
                                endpos = parse(Int,regexposstr[2])
                                index = 0
                                for posindex=startpos:endpos
                                    if postype == "N:"
                                        push!(alignment.noncoding, posindex)
                                    elseif postype == "1:" && (!modulo3 || index % 3 == 0)
                                        push!(alignment.coding1, posindex)
                                    elseif postype == "2:"  && (!modulo3 || index % 3 == 0)
                                        push!(alignment.coding2, posindex)
                                    elseif postype == "3:"  && (!modulo3 || index % 3 == 0)
                                        push!(alignment.coding3, posindex)
                                    end
                                    index += 1
                                end
                                #println(regexposstr)
                            elseif length(strip(posstr)) > 0
                                if postype == "N:"
                                    push!(alignment.noncoding, parse(Int,posstr))
                                elseif postype == "1:"
                                    push!(alignment.coding1, parse(Int,posstr))
                                elseif postype == "2:"
                                    push!(alignment.coding2, parse(Int,posstr))
                                elseif postype == "3:"
                                    push!(alignment.coding3, parse(Int,posstr))
                                end
                            end
                        end
                        #println(postype)
                        #println(posspl)
                    end
                    #println(regexcodon)
                end
            end

        end
    end
    datalen = length(alignment.sequences[1])
    if length(alignment.noncoding) == 0 && length(alignment.coding1) == 0 && length(alignment.coding2) == 0 && length(alignment.coding3) == 0
        for i=1:datalen
            push!(alignment.noncoding, i)
        end
    end
    alignment.codonindices = getcodonindices(datalen, alignment.noncoding, alignment.coding1, alignment.coding2, alignment.coding3)
    close(fin)
    return alignment
end

function getcodonindices(len::Int, noncoding::Array{Int,1}, coding1::Array{Int,1}, coding2::Array{Int,1}, coding3::Array{Int,1})
    codonindices = zeros(Int,len)
    codonindex = 1
    for (c1,c2,c3) in zip(coding1,coding2,coding3)
        codonindices[c1] = codonindex
        codonindices[c2] = codonindex
        codonindices[c3] = codonindex
        codonindex += 1
    end
    return codonindices
end
#=

alignment = Alignment()
push!(alignment.sequences, "AAA")
push!(alignment.sequences, "AAT")
push!(alignment.taxa, "cat")
push!(alignment.taxa, "dog")
alignment.charactertype = 1
push!(alignment.noncoding, 1)
push!(alignment.noncoding, 2)
push!(alignment.noncoding, 3)
push!(alignment.noncoding, 4)
push!(alignment.coding1, 5)
push!(alignment.coding2, 6)
push!(alignment.coding3, 7)
push!(alignment.coding1, 8)
push!(alignment.coding2, 9)
push!(alignment.coding3, 10)
push!(alignment.treestrings, "(cat,dog);")
#savenexus(alignment, "test.nex")
println(loadnexus("test.nex"))=#
#println(loadnexus("test.nex"))

#=
println(loadnexus("lysin.nex"))=
savenexus(loadnexus("lysin.nex"), "test2.nex")
=#
#println(loadnexus("lysin.nex"))
