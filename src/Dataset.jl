using SparseArrays

include("AlignmentIO.jl")
include("Tree.jl")

nucmapping = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, 'U' => 4)

nongaps = "ACGTUBDHKMRSVWY"
function isgap(char::Char)
    return !(char in nongaps)
end

function setnucleotides(data::Array{Float64,3}, seq::Int, col::Int, char::Char)
  c = uppercase(char)
  if c == 'A'
    data[seq,col,1] = 1.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 0.0
  elseif c == 'C'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 0.0
  elseif c == 'G'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 0.0
  elseif c == 'T'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 1.0
  elseif c == 'U'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 1.0
  elseif c == 'B'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 1.0
  elseif c == 'D'
    data[seq,col,1] = 1.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 1.0
  elseif c == 'H'
    data[seq,col,1] = 1.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 1.0
  elseif c == 'K'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 1.0
  elseif c == 'M'
    data[seq,col,1] = 1.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 0.0
  elseif c == 'R' # purine
    data[seq,col,1] = 1.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 0.0
  elseif c == 'S'
    data[seq,col,1] = 0.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 0.0
  elseif c == 'V'
    data[seq,col,1] = 1.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 0.0
  elseif c == 'W'
    data[seq,col,1] = 1.0
    data[seq,col,2] = 0.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 1.0
  elseif c == 'Y' # pyrimidine
    data[seq,col,1] = 0.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 0.0
    data[seq,col,4] = 1.0
  else
    data[seq,col,1] = 1.0
    data[seq,col,2] = 1.0
    data[seq,col,3] = 1.0
    data[seq,col,4] = 1.0
  end
end

export Dataset
mutable struct Dataset
  data::Array{Float64,3} # sequence by column by alphabet
  logdata::Array{Float64,3}  # sequence by column by alphabet
  datalist::Array{Array{Float64,2},1}
  logdatalist::Array{Array{Float64,2},1}
  obsfreqs::Array{Float64,1}
  numseqs::Int
  numcols::Int
  numnodes::Int
  seqnametoindex::Dict{AbstractString,Int}
  root::TreeNode
  subcolumnrefs::Array{Int,2}
  gapfrequency::Array{Float64,1}
  sequences::Array{AbstractString,1}

  noncoding::Array{Int,1}
  coding1::Array{Int,1}
  coding2::Array{Int,1}
  coding3::Array{Int,1}
  codonindices::Array{Int,1}
  maxbasepairprobs::SparseMatrixCSC{Float64,Int64}

  function Dataset(alignment::Alignment, doroottree::Bool=true)
    sequences = alignment.sequences
    names = alignment.taxa
    seqnametoindex = Dict{AbstractString,Int}()

    len = 0
    seqindex = 1
    for (taxon,sequence) in zip(alignment.taxa, alignment.sequences)
        len = length(sequence)
        seqnametoindex[taxon] = seqindex
        seqindex += 1
    end

    numseqs = length(sequences)
    numcols = len
    data = zeros(Float64,numseqs,numcols,4)
    obsfreqs = zeros(Float64,4)
    gapfrequency = zeros(Float64,numcols)
    for s=1:numseqs
      seq = sequences[s]
      for j=1:numcols
        nuc = get(nucmapping,seq[j],0)
        if nuc > 0
          #data[s,j,nuc] = 1.0
          obsfreqs[nuc] += 1.0
        else
          #data[s,j,:] = 1.0
          gapfrequency[j] += 1.0
        end
        setnucleotides(data, s, j, seq[j])
      end
    end

    datalist = Array{Float64,2}[]
    logdatalist = Array{Float64,2}[]
    for j=1:numcols
      mat = zeros(Float64, numseqs, 4)
      for s=1:numseqs
        for a=1:4
          mat[s,a] = data[s,j,a]
        end
      end
      push!(datalist, mat)
      push!(logdatalist, log.(mat))
    end

    gapfrequency /= numseqs

    #=
    newickin = open(treefile,"r")
    newickstring = strip(readlines(newickin)[1])
    close(newickin)=#
    newickstring = alignment.treestrings[1]
    root = gettreefromnewick(newickstring)
    if doroottree
	     root = roottree(gettreefromnewick(newickstring), 1)
    end
    annotatetree(root,seqnametoindex)

    nodelist = getnodelist(root)
    numnodes = length(nodelist)
    subcolumnrefs = zeros(Int,numnodes,numcols)
    s = 0
    for node in nodelist
      dict = Dict()
      for col=1:numcols
        pattern = (node.nodeindex, getpattern(data, node, col))
        if haskey(dict, pattern)
          subcolumnrefs[node.nodeindex,col] = dict[pattern]
          s += length(pattern)-1
        else
          subcolumnrefs[node.nodeindex,col] = col
          dict[pattern] = col
        end
      end
    end
    return new(data, log.(data), datalist,logdatalist, obsfreqs/sum(obsfreqs), numseqs, numcols, numnodes, seqnametoindex, root,subcolumnrefs, gapfrequency,sequences, alignment.noncoding, alignment.coding1, alignment.coding2, alignment.coding3, alignment.codonindices, spzeros(numcols, numcols))
  end
end

export hascanonicalbasepair
function hascanonicalbasepair(dataset::Dataset, x::Int, y::Int)
    pairs = Tuple{Int,Int}[(2,3),(3,2),(1,4),(4,1),(3,4),(4,3)]

    for pair in pairs
        p1 = pair[1]
        p2 = pair[2]
        for s=1:dataset.numseqs
          if dataset.data[s,x,p1] > 0.0 && dataset.data[s,y,p2] > 0.0
              if !isgap(dataset.sequences[s][x]) && !isgap(dataset.sequences[s][y])
                  return true
              end
          end
        end
    end
    return false
end

export getindexbyseqcol
function getindexbyseqcol(dataset::Dataset, seq::Int, col::Int)
  return dataset.datalist[col][seq,:]
end

export fillindexbyseqcol
function fillindexbyseqcol(dataset::Dataset, seq::Int, col::Int, v::Array{Float64,1})
  maxindex = length(dataset.datalist[col][seq,:])
  len = max(length(v), maxindex)
  for a=1:len
    if a <= maxindex
      v[a] = dataset.datalist[col][seq,a]
    else
      v[a] = 0.0
    end
  end
end

function getdatasetfromalignment(fastafile::AbstractString, newickfile::AbstractString)
    alignment = getalignmentfromfastaandnewick(fastafile, newickfile)
    return Dataset(alignment, true)
end
