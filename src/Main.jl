push!(LOAD_PATH,string(@__DIR__))
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using LG
using CommonUtils
using FastaIO
using NLopt
using Formatting
using Distributions
using StatsBase
using ArgParse
using CTMCs
using JSON
using Viridis
using HypothesisTests
using ChiSquaredEvolution
using Random
using Printf
using Nullables
using LinearAlgebra

aminoacids = "ACDEFGHIKLMNPQRSTVWY"

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--alignment"
          help = "alignment file to be used in FASTA format"
          arg_type = String
          required = true
        "--tree"
          help = "tree file to be used in Newick format"
          arg_type = String
          required = false
        "--annotations"
          help = ""
          arg_type = String
          required = true
    end

    return parse_args(s)
end

function getTmatrix(tau::Float64, traitfreqs::Array{Float64,1})
    HP = 2
    LP = 1
    T = zeros(2,2)
    for m=1:2
        for n=1:2
            if m != n
                T[m,n] = traitfreqs[n]*tau
                T[m,m] -= T[m,n]
            end
        end
    end

    return T, traitfreqs
end

mutable struct AATraitParameters
    mu::Float64
    lambda::Float64
    tau::Float64
    p::Float64
    traitfreqs::Array{Float64,1}

    function AATraitParameters(mu::Float64, lambda::Float64, tau::Float64, p::Float64)
        new(mu, lambda, tau, p, Float64[p, 1.0-p])
    end

    function AATraitParameters(arr::Array{Float64,1})
        new(arr[1], arr[2], arr[3], arr[4], Float64[arr[4], 1.0-arr[4]])
    end
end

function getarray(params::AATraitParameters)
    return Float64[params.mu,params.lambda,params.tau,params.p]
end

function getAATraitmatrix(params::AATraitParameters, A::Array{Float64,2}, aafreqs::Array{Float64,1}, targetaa::Int)
    mu = params.mu
    lambda = params.lambda
    tau = params.tau
    traitfreqs = params.traitfreqs

    Alen =size(A,1)
    HP = 2
    LP = 1
    Q = zeros(2*Alen,2*Alen)
    for m=1:2
        for n=1:2
            for aa1=1:Alen
                for aa2=1:Alen
                    i = (m-1)*Alen + aa1
                    j = (n-1)*Alen + aa2
                    if i != j
                        if aa1 != aa2 && m == n
                            if n == HP
                                if aa2 == targetaa
                                    Q[i,j] = mu*A[aa1,aa2]*lambda
                                elseif aa2 != targetaa
                                    Q[i,j] = mu*A[aa1,aa2]/lambda
                                end
                            elseif n == LP
                                Q[i,j] = mu*A[aa1,aa2]
                            end
                        elseif m != n && aa1 == aa2
                            if n == HP
                                if aa2 == targetaa
                                    Q[i,j] = traitfreqs[n]*tau*lambda
                                elseif aa2 != targetaa
                                    Q[i,j] = traitfreqs[n]*tau/lambda
                                end
                            elseif n == LP
                                Q[i,j] = traitfreqs[n]*tau
                            end
                        end
                        Q[i,i] -= Q[i,j]
                    end
                end
            end
        end
    end
    freqs = zeros(Float64, 2*Alen)
    for n=1:2
        for aa=1:Alen
            i = (n-1)*Alen + aa
            if aa == targetaa && n == LP
                freqs[i] = aafreqs[aa]*traitfreqs[n]
            elseif aa != targetaa && n == LP
                freqs[i] = aafreqs[aa]*traitfreqs[n]
            elseif aa == targetaa && n == HP
                freqs[i] = aafreqs[aa]*traitfreqs[n]*lambda
            elseif aa != targetaa && n == HP
                freqs[i] = aafreqs[aa]*traitfreqs[n]/lambda
            end
        end
    end
    freqs /= sum(freqs)
    marginaltraitfreqs = zeros(Float64,2)
    for n=1:2
        for aa=1:Alen
            i = (n-1)*Alen + aa
            marginaltraitfreqs[n] += freqs[i]
        end
    end
    #println("A",sum(diagm(freqs)*Q,1))
    v = freqs'*Q
    for i=1:40
        if abs(v[i]) >= 1e-10
            println("lambda", lambda)
            println("not stationary: ", i,",",abs(v[i]))
            break
        end
    end
    for i=1:40
        for j=1:40
            delta = freqs[i]*Q[i,j] - freqs[j]*Q[j,i]
            if abs(delta) >= 1e-10
                println("lambda", lambda)
                println("not reversible: ", i,",",j,",",delta)
                break
            end
        end
    end

    return Q, freqs, marginaltraitfreqs
end

function getaacolumn(sequences::Array{AbstractString,1}, col::Int)
    ret = zeros(Int, length(sequences))
    s = 1
    for seq in sequences
        aa = string(uppercase(seq[col]))
        index = indexof(aa,aminoacids)
        ret[s] = index
        s += 1
    end
    return ret
end

function felsensteinstack(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, logm::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, alphabet::Int)

    for node in nodelist
        if !isleafnode(node)
            for a=1:alphabet
                likelihoods[node.nodeindex,a] = -Inf
            end
        end
    end
  stack = Int[1]
  while length(stack) > 0
    nodeindex = stack[end]
    node = nodelist[nodeindex]
    if isleafnode(node)
        pop!(stack)
    else
        leftchildindex = node.children[1].nodeindex
        rightchildindex = node.children[2].nodeindex

        cont = true
        if likelihoods[leftchildindex, 1] == -Inf
          push!(stack, leftchildindex)
          cont = false
        end
        if likelihoods[rightchildindex, 1] == -Inf
          push!(stack, rightchildindex)
          cont = false
        end

        if cont
          v = (transprobs[leftchildindex]*likelihoods[leftchildindex,:]).*(transprobs[rightchildindex]*likelihoods[rightchildindex,:])
          likelihoods[nodeindex,:] = (transprobs[leftchildindex]*likelihoods[leftchildindex,:]).*(transprobs[rightchildindex]*likelihoods[rightchildindex,:])
          m = maximum(likelihoods[nodeindex,:])
          if m < 1e-20
            likelihoods[nodeindex,:] /= m
            logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
          else
            logm[nodeindex] = logm[leftchildindex] + logm[rightchildindex]
          end
          pop!(stack)
        end
    end
  end

  return likelihoods
end

function backwardssampling(rng::AbstractRNG, nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, Q::Array{Float64,2}, freqs::Array{Float64,1}, nsamples::Int=1)
    logm = zeros(Float64, size(likelihoods,1))
    transprobs = gettransprobmatrices(nodelist, Q)
    felsensteinstack(nodelist, likelihoods,logm, transprobs, size(Q,1))
    rootliks = likelihoods[1,:].*freqs
    rootliks /= sum(rootliks)
    samples = zeros(Int, nsamples, length(nodelist))
    for s=1:nsamples
        samples[s, 1] = CommonUtils.sample(rng, rootliks)
        stack = Int[1]
        while length(stack) > 0
          nodeindex = stack[end]
          node = nodelist[nodeindex]
          parentc = samples[s, nodeindex]
          pop!(stack)
          for child in node.children
              if samples[s, child.nodeindex] == 0
                  samples[s, child.nodeindex] = CommonUtils.sample(rng, transprobs[child.nodeindex][parentc,:].*likelihoods[child.nodeindex,:])
                  if !isleafnode(child)
                      push!(stack, child.nodeindex)
                  end
              end
          end
      end
    end
    return samples
end

function gettransprobmatrices(nodelist::Array{TreeNode,1}, Q::Array{Float64,2})
    transprobs = Array{Float64,2}[]
    decomposition = eigen(Q)
    D, V = decomposition.values, decomposition.vectors
    Vi = inv(V)
    for node in nodelist
        Pt = V*Diagonal(exp.(D*node.branchlength))*Vi
        push!(transprobs, Pt)
    end
    return transprobs
end

function prettyprintmatrix(mat::Array{Float64,2})
    for i=1:size(mat,2)
        print(lpad(string(i), 6, '0'),"\t")
    end
    println()
    for i=1:size(mat,1)
        for j=1:size(mat,2)
            print(@sprintf("%0.4f", mat[i,j]), "\t")
        end
        println("")
    end
end

function computesiteloglikelihoodhelper(param::Array{Float64,1}, datalikelihoods::Array{Float64,2}, mles::Array{Float64,2}, mlparams::Dict{Tuple{Int,Int}, AATraitParameters}, nodelist::Array{TreeNode,1}, data::Array{Float64,3}, subcolumnrefs::Array{Int,2}, initialparams::AATraitParameters, targetaa::Int, col::Int, maxallsiteiters::Int)
    try
        global optiters
        optiters += 1
        if (optiters-1) % 50 == 0
            return computesiteloglikelihoods(mles, mlparams, nodelist, data, subcolumnrefs, AATraitParameters(param), targetaa)[col]
        else
            return getaatraitcolumnloglikelihood(nodelist, datalikelihoods, AATraitParameters(param), targetaa)
        end
    catch e
        return -1e20
    end
end

function optimizeaatraitmodel(datalikelihoods::Array{Float64,2},mles::Array{Float64,2}, mlparams::Dict{Tuple{Int,Int}, AATraitParameters}, nodelist::Array{TreeNode,1}, data::Array{Float64,3}, subcolumnrefs::Array{Int,2}, initialparams::AATraitParameters, targetaa::Int, col::Int, nloptmethod::Symbol=:LN_COBYLA, fixtraitparams::Bool=false)
     initialarr = getarray(initialparams)
     println("A",initialarr)
     global optiters
     optiters = 0
     maxoptiter = 2000
     opt = Opt(nloptmethod, 4)
    #localObjectiveFunction = ((param, grad) -> computesiteloglikelihoods(mles, mlparams, nodelist, data, subcolumnrefs, AATraitParameters(param), targetaa)[col])
    localObjectiveFunction = ((param, grad) -> computesiteloglikelihoodhelper(param, datalikelihoods, mles, mlparams, nodelist, data, subcolumnrefs, initialparams, targetaa, col, 100))
    #=
    if opt
        localObjectiveFunction = ((param, grad) -> getaatraitcolumnloglikelihood(nodelist, likelihoods, AATraitParameters(param), targetaa))
    end=#
    lower = ones(Float64, 4)*1e-3
    upper = ones(Float64, 4)*1e3
    if targetaa == 0
        lower[2] = 1.0
        upper[2] = 1.0
    end
    lower[4] = 0.01
    upper[4] = 0.99
    if fixtraitparams
        lower[3] = initialarr[3]
        upper[3] = initialarr[3]
        lower[4] = initialarr[4]
        upper[4] = initialarr[4]
    end
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    if mode == 1
        xtol_rel!(opt,1e-3)
        maxeval!(opt, maxoptiter)
    else
        xtol_rel!(opt,1e-6)
        maxeval!(opt, 3*maxoptiter)
    end
    max_objective!(opt, localObjectiveFunction)
    println("B",initialarr)
    initialarr[4] = max(0.011,initialarr[4])
    initialarr[4] = min(0.989,initialarr[4])
    (maxll,maxparams,ret) = optimize(opt, initialarr)
    println("Optimal: ", maxll,"\t",maxparams)
    return maxll,maxparams
end

function optimizetraitmodel(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, initialparams::AATraitParameters, nloptmethod::Symbol=:LN_COBYLA, maxoptiter::Int=20000)
    opt = Opt(nloptmethod, 2)
    localObjectiveFunction = ((param, grad) -> gettraitcolumnloglikelihood(nodelist, likelihoods, param[1], param[2]))
    lower = ones(Float64, 2)*1e-3
    upper = ones(Float64, 2)*1e3
    lower[2] = 0.01
    upper[2] = 0.99
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, maxoptiter)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, getarray(initialparams)[3:4])
    finalparams = getarray(initialparams)
    finalparams[3] = minx[1]
    finalparams[4] = minx[2]
    return AATraitParameters(finalparams)
end

function getaalikelihoods(nodelist::Array{TreeNode,1}, sequences::Array{AbstractString,1}, aacolumn::Array{Int,1}, seqindextonodeindex::Array{Int,1})
    alphabet = 20
    likelihoods = ones(length(nodelist),alphabet)*-Inf
    for seqindex=1:length(sequences)
        nodeindex = seqindextonodeindex[seqindex]
        for a=1:alphabet
            likelihoods[nodeindex,a] = 0.0
        end
        if aacolumn[seqindex] == 0
            for a=1:alphabet
                likelihoods[nodeindex,a] = 1.0
            end
        else
            likelihoods[nodeindex,aacolumn[seqindex]] = 1.0
        end
    end
    return likelihoods
end

function getaatraitlikelihoods(nodelist::Array{TreeNode,1}, sequences::Array{AbstractString,1}, aacolumn::Array{Int,1}, traitcolumn::Array{Int,1}, seqindextonodeindex::Array{Int,1})
    alphabet = 40
    likelihoods = ones(length(nodelist),alphabet)*-Inf
    for seqindex=1:length(sequences)
        nodeindex = seqindextonodeindex[seqindex]
        for a=1:alphabet
            likelihoods[nodeindex,a] = 0.0
        end
        if traitcolumn[seqindex] == 0 && aacolumn[seqindex] == 0
            for a=1:alphabet
                likelihoods[nodeindex,a] = 1.0
            end
        elseif traitcolumn[seqindex] == 0
            likelihoods[nodeindex,aacolumn[seqindex]] = 1.0
            likelihoods[nodeindex,20+aacolumn[seqindex]] = 1.0
        elseif  aacolumn[seqindex] == 0
            for aa=1:20
                likelihoods[nodeindex,(traitcolumn[seqindex]-1)*20+aa] = 1.0
            end
        else
            likelihoods[nodeindex,(traitcolumn[seqindex]-1)*20 + aacolumn[seqindex]] = 1.0
        end
    end
    return likelihoods
end

function gettraitlikelihoods(nodelist::Array{TreeNode,1}, traitcolumn::Array{Int,1}, seqindextonodeindex::Array{Int,1})
    likelihoods = ones(length(nodelist),2)*-Inf
    for seqindex=1:length(traitcolumn)
        nodeindex = seqindextonodeindex[seqindex]
        for a=1:size(likelihoods,2)
            likelihoods[nodeindex,a] = 0.0
        end
        if traitcolumn[seqindex] == 0
            for a=1:size(likelihoods,2)
                likelihoods[nodeindex,a] = 1.0
            end
        else
            likelihoods[nodeindex,traitcolumn[seqindex]] = 1.0
        end
    end
    return likelihoods
end

function getaatraitcolumnloglikelihood(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, params::AATraitParameters, targetaa::Int)
    alphabet = 40
    Q,freqs,marginaltraitfreqs = getAATraitmatrix(params, LGmatrix, LGfreqs, targetaa)
    logm = zeros(Float64, size(likelihoods,1))
    felsensteinstack(nodelist, likelihoods,logm, gettransprobmatrices(nodelist, Q), alphabet)
    totalloglikelihood = logm[1]+log.(dot(likelihoods[1,:], freqs))
    #println(totalloglikelihood,"\t",mu,"\t",lambda,"\t",tau,"\t",p,"\t",marginaltraitfreqs)
    #println(logm)
    return totalloglikelihood
end

function getaacolumnloglikelihood(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, mu::Float64)
    logm = zeros(Float64, size(likelihoods,1))
    felsensteinstack(nodelist, likelihoods,logm, gettransprobmatrices(nodelist, mu*LGmatrix), 20)
    totalloglikelihood = logm[1]+log.(dot(likelihoods[1,:], LGfreqs))
    return totalloglikelihood
end

function gettraitcolumnloglikelihood(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2},tau::Float64, p::Float64)
    T,freqs = getTmatrix(tau, Float64[p, 1.0-p])
    logm = zeros(Float64, size(likelihoods,1))
    felsensteinstack(nodelist, likelihoods,logm, gettransprobmatrices(nodelist, T), 2)
    totalloglikelihood = logm[1]+log.(dot(likelihoods[1,:], freqs))
    return totalloglikelihood
end

function getaatraitdata(nodelist::Array{TreeNode,1}, seqindextonodeindex::Array{Int,1}, sequences::Array{AbstractString,1}, traits::Array{Int,1})
    numseqs = length(sequences)
    numcols = length(sequences[1])
    alphabet = 40
    data = zeros(Float64,numseqs,numcols,alphabet)
    for seqindex=1:length(sequences)
        seq = sequences[seqindex]
        for col=1:numcols
            aa = max(0, indexof(string(uppercase(seq[col])),aminoacids))
            nodeindex = seqindextonodeindex[seqindex]
            if traits[seqindex] == 0 && aa == 0
                for a=1:alphabet
                    data[seqindex,col,a] = 1.0
                end
            elseif traits[seqindex] == 0
                data[seqindex,col,aa] = 1.0
                data[seqindex,col,20+aa] = 1.0
            elseif  aa == 0
                for aa=1:20
                    data[seqindex,col,(traits[seqindex]-1)*20+aa] = 1.0
                end
            else
                data[seqindex,col,(traits[seqindex]-1)*20 + aa] = 1.0
            end
        end
    end


    subcolumnrefs = zeros(Int,length(nodelist),numcols)
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
    return data,subcolumnrefs
end

function getaadata(nodelist::Array{TreeNode,1}, seqindextonodeindex::Array{Int,1}, sequences::Array{AbstractString,1})
    numseqs = length(sequences)
    numcols = length(sequences[1])
    alphabet = 20
    data = zeros(Float64,numseqs,numcols,alphabet)
    for seqindex=1:length(sequences)
        seq = sequences[seqindex]
        for col=1:numcols
            aa = max(0, indexof(string(uppercase(seq[col])),aminoacids))
            nodeindex = seqindextonodeindex[seqindex]
            if aa == 0
                for a=1:alphabet
                    data[seqindex,col,a] = 1.0
                end
            else
                data[seqindex,col,aa] = 1.0
            end
        end
    end


    subcolumnrefs = zeros(Int,length(nodelist),numcols)
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
    return data,subcolumnrefs
end

function felsensteinstack(nodelist::Array{TreeNode,1}, data::Array{Float64,3}, subcolumnrefs::Array{Int,2}, col::Int, likelihoods::Array{Float64,3}, logm::Array{Float64,2}, transprobs::Array{Array{Float64,2},1}, alphabet::Int)
  stack = Int[1]
  while length(stack) > 0
    nodeindex = stack[end]
    node = nodelist[nodeindex]

    subcol = subcolumnrefs[nodeindex,col]
    if subcol < col
        for a=1:alphabet
            likelihoods[nodeindex,col,a] = likelihoods[nodeindex,subcol,a]
        end
        logm[col,nodeindex] = logm[subcol,nodeindex]
        pop!(stack)
    elseif isleafnode(node)
        for a=1:alphabet
            likelihoods[nodeindex,col,a] = data[node.seqindex,col,a]
        end
        pop!(stack)
    else
        leftchildindex = node.children[1].nodeindex
        rightchildindex = node.children[2].nodeindex

        cont = true
        if likelihoods[leftchildindex,col, 1] == -Inf
          push!(stack, leftchildindex)
          cont = false
        end
        if likelihoods[rightchildindex,col, 1] == -Inf
          push!(stack, rightchildindex)
          cont = false
        end

        if cont
          likelihoods[nodeindex,col,:] = (transprobs[leftchildindex]*likelihoods[leftchildindex,col,:]).*(transprobs[rightchildindex]*likelihoods[rightchildindex,col,:])
          m = maximum(likelihoods[nodeindex,col,:])
          if m < 1e-20
            likelihoods[nodeindex,col,:] /= m
            logm[col,nodeindex] = log(m) + logm[col,leftchildindex] + logm[col,rightchildindex]
          else
            logm[col,nodeindex] = logm[col,leftchildindex] + logm[col,rightchildindex]
          end
          pop!(stack)
        end
    end
  end

  return likelihoods
end

optiters = 0
function computesiteloglikelihoods(mles::Array{Float64,2}, mlparams::Dict{Tuple{Int,Int}, AATraitParameters}, nodelist::Array{TreeNode,1}, data::Array{Float64,3}, subcolumnrefs::Array{Int,2}, params::AATraitParameters, targetaa::Int)
    tic()
    numcols = size(data,2)
    siteloglikelihoods = zeros(Float64, numcols)
    try
        likelihoods = ones(Float64,length(nodelist),numcols,40)*-Inf
        Q, freqs, marginaltraitfreqs = getAATraitmatrix(params, LGmatrix, LGfreqs, targetaa)
        transprobs = gettransprobmatrices(nodelist, Q)
        logm = zeros(Float64,numcols,length(nodelist))
        for col=1:numcols
            felsensteinstack(nodelist, data, subcolumnrefs, col, likelihoods, logm, transprobs, 40)
            if subcolumnrefs[1,col] < col
                siteloglikelihoods[col] = siteloglikelihoods[subcolumnrefs[1,col]]
            else
                siteloglikelihoods[col] = logm[col,1] + log.(dot(likelihoods[1,col,:],freqs))
            end
            t = targetaa
            if t == 0
                t = 21
                for s=1:21
                    if siteloglikelihoods[col] > mles[col,s]
                        mles[col,s] = siteloglikelihoods[col]
                        if s == 21
                            mlparams[(col,0)] = params
                        else
                            mlparams[(col,s)] = params
                        end
                    end
                end
            elseif siteloglikelihoods[col] > mles[col,t]
                mles[col,t] = siteloglikelihoods[col]
                mlparams[(col,targetaa)] = params
            end
        end
    catch e
        println(e)
        println(catch_stacktrace())
        #exit()
    end
    println(optiters)
    toc()
    return siteloglikelihoods
end

function performchi2association(outfile::AbstractString, sequences::Array{AbstractString,1}, traits::Array{Int,1}, annotationnames)
    numcols = length(sequences[1])
    numseqs = length(sequences)
    fout = open(outfile, "w")
    println(fout,"\"Position\",\"Amino acid\",\"Chi2\",\"Fisher's pvalue\",\"AA count\",\"notAA count\",\"$(annotationnames[1])+AA\",\"Exp[$(annotationnames[1])+AA]\",\"$(annotationnames[2])+AA\",\"Exp[$(annotationnames[2])+AA]\",\"$(annotationnames[1])-AA\",\"Exp[$(annotationnames[1])-AA]\",\"$(annotationnames[2])-AA\",\"Exp[$(annotationnames[2])-AA]\"")
    for i=1:numcols
        for targetaa=1:20
            table = zeros(Float64,2,2)
            for s=1:numseqs
                aaindex = 0
                if aminoacids[targetaa] == sequences[s][i]
                    aaindex = 2
                elseif sequences[s][i] != '-'
                    aaindex = 1
                end

                if aaindex > 0 && traits[s] > 0
                    table[aaindex,traits[s]] += 1.0
                end
            end
            expected, chi2 = expectedchi(table)
            fisherspvalue = fisherstest(table)
            if isnan(chi2)
                chi2 = 0.0
            elseif table[2,2] > expected[2,2]
                chi2 = -chi2
            end

            println(fout,i,",",aminoacids[targetaa],",",chi2,",", fisherspvalue,",", table[2,1]+table[2,2],",", table[1,1]+table[1,2],",", table[2,1],",",@sprintf("%0.2f", expected[2,1]),",",table[2,2],",",@sprintf("%0.2f", expected[2,2]),",", table[1,1],",",@sprintf("%0.2f",expected[1,1]),",",table[1,2],",",@sprintf("%0.2f",expected[1,2]))
        end
    end
    close(fout)
end

function main()
    parsed_args = parse_commandline()
    fastafile = parsed_args["alignment"]
    treefile = parsed_args["tree"]
    rng = MersenneTwister(1)
    zup = 1

    prioritycolumns = Int[143]
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            push!(sequences, seq)
        end
    end
    numseqs = length(sequences)
    numcols = length(sequences[1])
    seqnametoindex = Dict{AbstractString,Int}()
    len = 0
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        len = length(sequence)
        seqnametoindex[taxon] = seqindex
        seqindex += 1
    end
    annotationnames = split(parsed_args["annotations"],",")

    traits = Int[]
    for name in names
        match = false
        for annotationindex=1:length(annotationnames)
            if occursin(annotationnames[annotationindex], name)
                push!(traits,annotationindex)
                match = true
                break
            end
        end
        if !match
            push!(traits,0)
        end
    end
    println(traits)

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
        if node.branchlength <= 1e-7
            node.branchlength = 1e-7
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
    seqindextonodeindex = zeros(Int,length(sequences))
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        seqindextonodeindex[seqindex] = seqnametonodeindex[taxon]
        nodelist[seqnametonodeindex[taxon]].seqindex = seqindex
        seqindex += 1
    end

    data,subcolumnrefs = getaatraitdata(nodelist,seqindextonodeindex,sequences, traits)


    performchi2association(string(fastafile,".chi2.csv"), sequences, traits, annotationnames)
    #=
    tic()
    mlesin = ones(Float64, numcols, 21)*-Inf
    mlparamsin = Dict{Tuple{Int,Int}, AATraitParameters}()
    siteloglikelihoods = computesiteloglikelihoods(mlesin, mlparamsin, nodelist, data, subcolumnrefs, AATraitParameters(1.0,1.0,1.0,0.5), 1)
    toc()
    println(siteloglikelihoods)
    println(mlesin)=#

    #exit()
    #=
    markednodes = zeros(Int,length(nodelist))

    for node in nodelist
        println(node.nodeindex)
        traitcounts = zeros(Int, 2)
        total = 0
        stack = Int[node.nodeindex]
        while length(stack) > 0
            currentnode = nodelist[stack[end]]
            if isleafnode(currentnode)
                trait = traits[seqnametoindex[currentnode.name]]
                if trait > 0
                    traitcounts[trait] += 1
                end
                total += 1
            else
                push!(stack, currentnode.children[1].nodeindex)
                push!(stack, currentnode.children[2].nodeindex)
            end
            pop!(stack)
        end
        if total == traitcounts[1]
            markednodes[node.nodeindex] = 1
        end
        if total == traitcounts[2]
            markednodes[node.nodeindex] = 2
        end
    end
    println(markednodes)=#
    markednodes = zeros(Int,length(nodelist))
    traitcounts = ones(Int, length(nodelist), 2)*-1
    stack = Int[1]
    while length(stack) > 0
        nodeindex = stack[end]
        currentnode = nodelist[nodeindex]
        if isleafnode(currentnode)
            trait = traits[seqnametoindex[currentnode.name]]
            for a=1:size(traitcounts,2)
                traitcounts[nodeindex,a] = 0
            end
            if trait > 0
                traitcounts[nodeindex, trait] = 1
            end
            pop!(stack)
        else
            cont = true
            for child in currentnode.children
                if traitcounts[child.nodeindex,1] == -1
                    push!(stack, child.nodeindex)
                    cont = false
                end
            end
            if cont
                for a=1:size(traitcounts,2)
                    traitcounts[nodeindex,a] = 0
                end
                for child in currentnode.children
                    traitcounts[nodeindex,:] += traitcounts[child.nodeindex,:]
                end
                pop!(stack)
            end
        end
    end
    for i=1:length(nodelist)
        nodeindex = nodelist[i].nodeindex
        parentnode = nodelist[i].parent
        parentindex = 0
        if !isnull(parentnode)
            parentindex = get(parentnode).nodeindex
        end
        if !(parentindex == 0 || maximum(traitcounts[parentindex,:]) == sum(traitcounts[parentindex,:])) && maximum(traitcounts[nodeindex,:]) == sum(traitcounts[nodeindex,:])
            markednodes[nodeindex] = argmax(traitcounts[nodeindex,:])
        end
    end


    columns = Int[]
    append!(columns, prioritycolumns)
    for c=1:numcols
        if !in(c,columns)
            push!(columns,c)
        end
    end
    println(columns)
    mles = ones(Float64,numcols,21)*-Inf
    optimisationstatus = zeros(Int,numcols,21)
    performcompleteoptimisation = false
    mlparams = Dict{Tuple{Int,Int}, AATraitParameters}()
    chi2lower = zeros(numcols,21,3)
    chi2medians = zeros(numcols,21,3)
    chi2upper = zeros(numcols,21,3)
    chi2lower_hponly = zeros(numcols,21,3)
    chi2medians_hponly = zeros(numcols,21,3)
    chi2upper_hponly = zeros(numcols,21,3)
    fishertable = zeros(numcols,21,3)
    fishertable_hponly = zeros(numcols,21,3)
    initialparams = Float64[1.0,1.0,1.0,0.5]
    aminoacidatleaf = zeros(Int,numcols,20)
    for col=1:numcols
         ret = getaacolumn(sequences, col)
         for s=1:size(ret,1)
             if ret[s] > 0
                 aminoacidatleaf[col,ret[s]] += 1
             end
         end
    end

    datalikelihoods = getaatraitlikelihoods(nodelist, sequences, getaacolumn(sequences, 143), traits, seqindextonodeindex)
    aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, 143), seqindextonodeindex)
    traitlikelihoods = gettraitlikelihoods(nodelist, traits, seqindextonodeindex)
    tll1 = gettraitcolumnloglikelihood(nodelist, traitlikelihoods, 6.040222915, 0.24263616)
    tll2 = gettraitcolumnloglikelihood(nodelist, traitlikelihoods, 7.456053451, 0.394240064)

    println("a", tll1)
    println("b", tll2)
    p2 = AATraitParameters(Float64[31.47925915,1.0,6.040222915,0.24263616])
    p1 = AATraitParameters(Float64[31.47925915,1.0,7.456053451, 0.394240064])
    cll1 = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, p2, 0)
    cll2 = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, p1, 0)
    aall = getaacolumnloglikelihood(nodelist,aalikelihoods, 31.47925915)
    println("e",aall)
    println("c",cll1,"\t", tll1+aall, "\t",tll1+aall-cll1)
    println("d",cll2,"\t", tll2+aall, "\t",tll2+aall-cll2)
    println(tll1-tll2)
    println(cll1-cll2)

    #initial = AATraitParameters(initialparams)
    initial = deepcopy(p1)

    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_DIRECT)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_CRS2_LM)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :LN_COBYLA)
    println(getarray(initial))
    initialparams = getarray(initial)
    for z=1:zup
        for col in columns
            datalikelihoods = getaatraitlikelihoods(nodelist, sequences, getaacolumn(sequences, col), traits, seqindextonodeindex)
            for targetaa=0:20
                if performcompleteoptimisation  || (z == 1 && (targetaa == 0 || aminoacidatleaf[col, targetaa] > 0)) || (z == 2 && (targetaa != 0 && aminoacidatleaf[col, targetaa] == 0))
                    label = string("independent\tcol: ", col)
                    if targetaa > 0
                        label = string("target: ", aminoacids[targetaa],"\tcol: ", col)
                    end
                    println(label)
                    nloptmethod = :LN_COBYLA
                    if targetaa == 0
                        nloptmethod = :LN_COBYLA
                    end
                    #maxll,maxparams = optimizeaatraitmodel(nodelist, likelihoods, targetaa, AATraitParameters(initialparams), mode)
                    initial = AATraitParameters(initialparams)
                    key = (col,targetaa)
                    if haskey(mlparams, key)
                        initial = mlparams[key]
                    end
                    #maxll,maxparams = optimizeaatraitmodel(datalikelihoods,mles, mlparams, nodelist, data, subcolumnrefs, initial, targetaa, col, mode, targetaa == 0)
                    maxll,maxparams = optimizeaatraitmodel(datalikelihoods,mles, mlparams, nodelist, data, subcolumnrefs, initial, targetaa, col, nloptmethod, targetaa == 0)
                    mlparams[key] = AATraitParameters(maxparams)
                    Q, freqs, marginaltraitfreqs = getAATraitmatrix(AATraitParameters(maxparams), LGmatrix, LGfreqs, targetaa)
                    println(marginaltraitfreqs)
                    if targetaa == 0
                        mles[col,21] = maxll
                        optimisationstatus[col,21] = 1
                        initialparams[2:4] = maxparams[2:4]
                        samples = backwardssampling(rng,nodelist, datalikelihoods, Q, freqs, 100)
                        for t=1:20
                            chi2medians[col,t,1], chi2lower[col,t,1], chi2upper[col,t,1], chi2medians[col,t,2], chi2lower[col,t,2], chi2upper[col,t,2],  chi2medians[col,t,3], chi2lower[col,t,3],  chi2upper[col,t,3], fishertable[col,t,1],fishertable[col,t,2],fishertable[col,t,3] = chisquaredtest(nodelist, markednodes, traits, samples, t, true)
                            chi2medians_hponly[col,t,1], chi2lower_hponly[col,t,1], chi2upper_hponly[col,t,1], chi2medians_hponly[col,t,2], chi2lower_hponly[col,t,2], chi2upper_hponly[col,t,2],  chi2medians_hponly[col,t,3], chi2lower_hponly[col,t,3], chi2upper_hponly[col,t,3], fishertable_hponly[col,t,1], fishertable_hponly[col,t,2], fishertable_hponly[col,t,3]  = chisquaredtest(nodelist, markednodes, traits, samples, t, false)
                        end
                        for t=1:20
                            if mles[col,21] > mles[col,t]
                                mles[col,t] = mles[col,21]
                                mlparams[(col,t)] = mlparams[key]
                            end
                        end
                    else
                        optimisationstatus[col,targetaa] = 1
                        mles[col,targetaa] = maxll
                        D = 2.0*(mles[col,targetaa]-mles[col,21])
                        println("chi2=",@sprintf("%0.2f", D),"\tpval=",ccdf(Chisq(1), D))
                    end
                    if targetaa == 0
                        initialparams = copy(maxparams)
                    end
                end
            end

            outfile = open(string(fastafile,".results.update.csv"),"w")
            println(outfile,"\"Site\",\"Target AA\",\"Count at leaf\",\"mu\",\"tau\",\"p\",\"lambda\",\"Optimisation\",\"MLL\",\"AIC\",\"chi2\",\"pvalue\",\"X->notX vs X->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"X->notX vs X->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\",\"Fisher's pvalue\",\"lower\",\"upper\"")
            for c=1:numcols
                if optimisationstatus[c,21] == 1
                    for t=0:20
                        key = (c,t)
                        if haskey(mlparams, key)
                            params = mlparams[key]
                            if t == 0
                                optstatuslabel = "Partial"
                                if optimisationstatus[c,21] == 1
                                    optstatuslabel = "Complete"
                                end
                                D = 0.0
                                aic = 2.0*3.0 - 2.0*mles[c,21]
                                println(outfile,c,",", "independent,","-,", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel, ",",mles[c,21], ",",aic, ",", "0.0", ",", "0.0", ",0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0")
                            else
                                optstatuslabel = "Partial"
                                if optimisationstatus[c,t] == 1
                                    optstatuslabel = "Complete"
                                end
                                D = 2.0*(mles[c,t]-mles[c,21])
                                pval = ccdf(Chisq(1), D)
                                if params.lambda < 1.0
                                    D = -D
                                end
                                aic = 2.0*4.0 - 2.0*mles[c,t]
                                println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", fishertable_hponly[col,t,1], ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", fishertable_hponly[col,t,2], ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", fishertable_hponly[col,t,3], ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", fishertable[col,t,1], ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", fishertable[col,t,2], ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", fishertable[col,t,3], ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                            end
                        end
                    end
                end
            end
            close(outfile)
        end
    end
end

function logprior(params::AATraitParameters)
    return logpdf(Exponential(1.0),params.mu) + logpdf(Exponential(1.0),params.tau) + logpdf(Normal(0.0, 1.0), log(params.lambda))
end

function mcmc()
    parsed_args = parse_commandline()
    fastafile = parsed_args["alignment"]
    treefile = parsed_args["tree"]
    rng = MersenneTwister(1)
    zup = 1

    prioritycolumns = Int[143]
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            push!(sequences, seq)
        end
    end
    numseqs = length(sequences)
    numcols = length(sequences[1])
    seqnametoindex = Dict{AbstractString,Int}()
    len = 0
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        len = length(sequence)
        seqnametoindex[taxon] = seqindex
        seqindex += 1
    end
    annotationnames = split(parsed_args["annotations"],",")

    traits = Int[]
    for name in names
        match = false
        for annotationindex=1:length(annotationnames)
            if occursin(annotationnames[annotationindex], name)
                push!(traits,annotationindex)
                match = true
                break
            end
        end
        if !match
            push!(traits,0)
        end
    end
    println(traits)

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
        if node.branchlength <= 1e-7
            node.branchlength = 1e-7
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
    seqindextonodeindex = zeros(Int,length(sequences))
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        seqindextonodeindex[seqindex] = seqnametonodeindex[taxon]
        nodelist[seqnametonodeindex[taxon]].seqindex = seqindex
        seqindex += 1
    end

    data,subcolumnrefs = getaatraitdata(nodelist,seqindextonodeindex,sequences, traits)

    current = AATraitParameters(1.0, 1.0, 1.0, 0.5)
    proposed = deepcopy(current)

    targetaa = 17
    col = 148

    datalikelihoods = getaatraitlikelihoods(nodelist, sequences, getaacolumn(sequences, col), traits, seqindextonodeindex)



    rng = MersenneTwister(976100981111394)
    currentll = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, current, targetaa) + logprior(current)
    proposedll = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, proposed, targetaa) + logprior(proposed)
    moveweights = ones(Float64,4)/4.0

    fout = open("mcmc.log", "w")
    println(fout, "iter\tll\tmu\ttau\tlambda\tloglambda\tp")
    for z=1:100000
        m = CommonUtils.sample(rng, moveweights)
        valid = true
        if m == 1
            proposed.mu += randn(rng)*2.5
            if proposed.mu <= 0.0
                valid = false
            end
        elseif m  == 2
            proposed.lambda += randn(rng)*2.5
            if proposed.lambda <= 0.0
                valid = false
            end
        elseif m == 3
            proposed.tau += randn(rng)*2.5
            if proposed.tau <= 0.0
                valid = false
            end
        elseif m == 4
            proposed.p += randn(rng)*0.1
            proposed.traitfreqs[1] = proposed.p
            proposed.traitfreqs[2] = 1.0 - proposed.p
            if proposed.p <= 0.0 || proposed.p >= 1.0
                valid = false
            end
        end

        if valid
            datalikelihoods = getaatraitlikelihoods(nodelist, sequences, getaacolumn(sequences, col), traits, seqindextonodeindex)
            proposedll = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, proposed, targetaa)  + logprior(proposed)
            delta = proposedll-currentll
        end

        if valid && exp(delta) > rand(rng)
            currentll = proposedll
            current = deepcopy(proposed)
        else
            proposed = deepcopy(current)
        end

        if z % 25 == 0
            println(fout, "$(z-1)\t$(currentll)\t$(current.mu)\t$(current.tau)\t$(current.lambda)\t$(log(current.lambda))\t$(current.p)")
            flush(fout)
        end

    end

end

function discretizegamma(shape::Float64, scale::Float64, numcategories::Int)
  if numcategories == 1
    return Float64[1.0]
  else
    catwidth = 1.0 / numcategories
    vs = Float64[catwidth/2.0 + (i-1.0)*catwidth for i=1:numcategories]
    gammadist = Gamma(shape, scale)
    return Float64[quantile(gammadist, v) for v in vs]
  end
end

function aa_alignmentloglikelihood(nodelist::Array{TreeNode,1}, numcols::Int, data::Array{Float64,3},subcolumnrefs::Array{Int,2}, mu::Float64, shape::Float64, numcats::Int=10)
    try
        rates = discretizegamma(shape, 1.0/shape, numcats)
        println(rates)
        siteloglikelihoodsconditionals = ones(Float64, numcols, numcats)*-Inf
        siteloglikelihoods = ones(Float64, numcols)*-Inf
        cat = 1
        for rate in rates
            logprob = log(1.0/numcats)
            freqs = LGfreqs
            Q = mu*LGmatrix*rate
            likelihoods = ones(Float64, length(nodelist), numcols, 20)*-Inf
            transprobs = gettransprobmatrices(nodelist, Q)
            logm = zeros(Float64, numcols, length(nodelist))
            for col=1:numcols
                felsensteinstack(nodelist, data, subcolumnrefs, col, likelihoods, logm, transprobs, 20)
                ll = logm[col,1] + log.(dot(likelihoods[1,col,:],freqs))
                siteloglikelihoods[col] = logsumexp(logprob+siteloglikelihoods[col], ll)
                siteloglikelihoodsconditionals[col,cat] = logprob+siteloglikelihoods[col]
            end
            cat += 1
        end
        println(mu,"\t",shape,"\t",sum(siteloglikelihoods))
        for col=1:numcols
            siteloglikelihoodsconditionals[col,:] = exp.(siteloglikelihoodsconditionals[col,:].-siteloglikelihoods[col])
        end
        return sum(siteloglikelihoods), siteloglikelihoodsconditionals
    catch e
        println("ERROR ", mu, "\t", shape)
        return -Inf,  ones(Float64, numcols, numcats)*-Inf
    end
end

function optimizealignmentlikelihood(nodelist::Array{TreeNode,1}, numcols::Int, data::Array{Float64,3},subcolumnrefs::Array{Int,2}, mu::Float64, shape::Float64, numcats::Int, maxoptiter::Int)
    opt = Opt(:LN_COBYLA, 2)
    localObjectiveFunction = ((param, grad) -> aa_alignmentloglikelihood(nodelist, numcols, data, subcolumnrefs, exp(param[1]), exp(param[2]), numcats)[1])
    #lower = ones(Float64, 2)*1e-3
    #upper = ones(Float64, 2)*1e3
    lower = ones(Float64,2)*-7.0
    upper = ones(Float64,2)*7.0
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, maxoptiter)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, log.(Float64[mu,shape]))
    return minf, exp.(minx)
end

#=
function pathreconstruction()
    rng = MersenneTwister(2184104820249809)
    parsed_args = parse_commandline()
    fastafile = parsed_args["alignment"]
    treefile = parsed_args["tree"]
    zup = 1

    #fastafile = "example.fas"1
    #treefile = "example.nwk"
    #fastafile = "H7_full_protein.fas"
    #treefile = "H7_full_protein.nwk"
    #prioritycolumns = Int[1,2,3,4,143,136,274,384,438]
    #prioritycolumns = Int[1,2,3,453]
    prioritycolumns = Int[]
    #prioritycolumns = Int[25,26,31,135]
    #fastafile = "H7_Genome/PB1/H7_PB1_alnP.fasta"
    #treefile = "H7_Genome/PB1/H7_PB1_alnP.fasta.nwk"
    #prioritycolumns = Int[1,2,3,4,154,152]
    #fastafile = "H7_alnP.fasta.norm.fas"
    #treefile = "H7_alnP.fasta.norm.fas.nwk"
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            push!(sequences, seq)
        end
    end
    numseqs = length(sequences)
    numcols = length(sequences[1])
    seqnametoindex = Dict{AbstractString,Int}()
    len = 0
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        len = length(sequence)
        seqnametoindex[taxon] = seqindex
        seqindex += 1
    end

    traits = Int[]
    for name in names
        if contains(name, ".LP.") || contains(name, "/LP/")  || contains(name, "|LP|")
            push!(traits,1)
        elseif contains(name, ".HP.")  || contains(name, "/HP/")  || contains(name, "|HP|")
            push!(traits,2)
        else
            push!(traits,0)
        end
    end

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
    seqindextonodeindex = zeros(Int,length(sequences))
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        seqindextonodeindex[seqindex] = seqnametonodeindex[taxon]
        nodelist[seqnametonodeindex[taxon]].seqindex = seqindex
        seqindex += 1
    end

    println(getjsontree(nodelist, traits))
    #println(prettyprintstring(nodelist[1]))
    #exit()
    jsondict = Dict{Any,Any}()
    numsamples = 25

    aadata,subcolumnrefs = getaadata(nodelist,seqindextonodeindex,sequences)
    minf, minx = optimizealignmentlikelihood(nodelist, numcols, aadata, subcolumnrefs, 0.22184120634178817, 0.11740585348963008)
    maxll, siteloglikelihoodsconditionals = aa_alignmentloglikelihood(nodelist, numcols, aadata, subcolumnrefs, minx[1], minx[2])
    rates = discretizegamma(minx[2], 1.0/minx[2], 10)
    selcols = Int[143, 384, 274, 438, 136,402,335,287,152,198,214]
    for col in selcols
        jsondict[string("Column ", col)] = Dict{Any,Any}()
        jsondict[string("Column ", col)]["samples"] = Dict{Any,Any}()
        jsondict[string("Column ", col)]["alphabet"] = aminoacids
        for s=1:numsamples
            index = CommonUtils.sample(rng, siteloglikelihoodsconditionals[col,:])
            #println(col,"\t",index,"\t",rates[index])
            logm = zeros(Float64, length(nodelist))
            aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)
            Q = minx[1]*LGmatrix*rates[index]
            sample = backwardssampling(rng, nodelist, aalikelihoods, Q, LGfreqs, 1)[1,:]
            jsondict[string("Column ", col)]["samples"][s] = Dict{Any,Any}()
            jsondict[string("Column ", col)]["samples"][s]["paths"] = []
            jsondict[string("Column ", col)]["samples"][s]["times"] = []
            push!(jsondict[string("Column ", col)]["samples"][s]["paths"], Int[sample[1]])
            push!(jsondict[string("Column ", col)]["samples"][s]["times"], Float64[0.0])
            for node in nodelist
                if node.nodeindex != 1
                    parent = get(node.parent)
                    paths,times = modifiedrejectionsampling(rng, Q*node.branchlength, sample[parent.nodeindex], sample[node.nodeindex], nothing)
                    push!(jsondict[string("Column ", col)]["samples"][s]["paths"], paths)
                    push!(jsondict[string("Column ", col)]["samples"][s]["times"], Float64[parse(Float64, Formatting.fmt("0.3f", v)) for v in times])
                end
            end
        end
    end

    initial = AATraitParameters(Float64[31.47925915,1.0,6.040222915,0.24263616])
    traitlikelihoods = gettraitlikelihoods(nodelist, traits, seqindextonodeindex)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_DIRECT, 2000)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_CRS2_LM, 2000)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :LN_COBYLA, 2000)
    println(getarray(initial))
    traitlikelihoods = gettraitlikelihoods(nodelist, traits, seqindextonodeindex)
    tau = initial.tau
    p = initial.p
    T,freqs = getTmatrix(tau, Float64[p, 1.0-p])
    jsondict["Trait"] = Dict{Any,Any}()
    jsondict["Trait"]["samples"] = Dict{Any,Any}()
    jsondict["Trait"]["alphabet"] = "LH"
    gettraitcolumnloglikelihood(nodelist, traitlikelihoods, tau, p)
    traitsamples = backwardssampling(rng, nodelist, traitlikelihoods, T, freqs, 100)
    for s=1:size(traitsamples,1)
        sample = traitsamples[s,:]
        jsondict["Trait"]["samples"][s] = Dict{Any,Any}()
        jsondict["Trait"]["samples"][s]["paths"] = []
        jsondict["Trait"]["samples"][s]["times"] = []
        push!(jsondict["Trait"]["samples"][s]["paths"], Int[sample[1]])
        push!(jsondict["Trait"]["samples"][s]["times"], Float64[0.0])
        for node in nodelist
            if node.nodeindex != 1
                parent = get(node.parent)
                paths,times = modifiedrejectionsampling(rng, T*node.branchlength, sample[parent.nodeindex], sample[node.nodeindex], nothing)
                push!(jsondict["Trait"]["samples"][s]["paths"], paths)
                push!(jsondict["Trait"]["samples"][s]["times"], Float64[parse(Float64, Formatting.fmt("0.3f", v)) for v in times])
            end
        end
    end
    columnorder = ["Trait"]
    for col in selcols
        push!(columnorder, string("Column ", col))
    end
    jsondict["columnorder"] = columnorder
    #jsondict["Trait"]["samples"][s]["times"] = AbstractString[Formatting.fmt(".4f", v) for v in jsondict["Trait"]["samples"][s]["times"]]

    ts = linspace(0.0, 1.0, 10)

    gradient = viridisgradient(_viridis_data)
    for node in nodelist
        nodescore = zeros(Float64, length(ts))
        for s=1:numsamples
            i = 1
            for t in ts
                sequence = Int[]
                for selcol in selcols
                    path = jsondict[string("Column ", selcol)]["samples"][s]["paths"][node.nodeindex]
                    time = jsondict[string("Column ", selcol)]["samples"][s]["times"][node.nodeindex]
                    push!(sequence, getchar(path, time, t))
                end
                nodescore[i] += HAscore(sequence)
                i += 1
            end
        end
        nodescore /= numsamples
        text = ""
        text = string(text, """<linearGradient id="branchgradient$(node.nodeindex)" x1="0%" y1="0%" x2="100%" y2="0%">\n""")
        for i=1:length(ts)
            rgba = getcolor(gradient, nodescore[i])
            text = string(text, """<stop offset="$(trunc(Int, ts[i]*100.0))%" style="stop-color:rgb($(rgba[1]),$(rgba[2]),$(rgba[3]));stop-opacity:1" />\n""")
        end
        text = string(text, "</linearGradient>\n")

        text = string(text, """<linearGradient id="initialgradient$(node.nodeindex)" x1="0%" y1="0%" x2="100%" y2="0%">\n""")
        rgba = getcolor(gradient, nodescore[end])
        text = string(text, """<stop offset="0%" style="stop-color:rgb($(rgba[1]),$(rgba[2]),$(rgba[3]));stop-opacity:1" />\n""")
        text = string(text, """<stop offset="100%" style="stop-color:rgb($(rgba[1]),$(rgba[2]),$(rgba[3]));stop-opacity:1" />\n""")
        text = string(text, "</linearGradient>\n")
        println(text)
    end


    println(JSON.json(jsondict))
end=#
function bhattacharya()
    parsed_args = parse_commandline()
    return bhattacharya(parsed_args["alignment"], parsed_args["tree"], split(parsed_args["annotations"],","), "outfile.txt")
end

function bhattacharya(fastafile::AbstractString, treefile::AbstractString, annotationnames, outfile::AbstractString)
    println("bhattacharya: ", fastafile)
    rng = MersenneTwister(2184104820249809)
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            push!(sequences, seq)
        end
    end
    numseqs = length(sequences)
    numcols = length(sequences[1])
    seqnametoindex = Dict{AbstractString,Int}()
    len = 0
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        len = length(sequence)
        seqnametoindex[taxon] = seqindex
        seqindex += 1
    end

    traits = Int[]
    for name in names
        match = false
        for annotationindex=1:length(annotationnames)
            if occursin(annotationnames[annotationindex], name)
                push!(traits,annotationindex)
                match = true
                break
            end
        end
        if !match
            push!(traits,0)
        end
    end
    println(traits)

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
        if node.branchlength <= 1e-7
            node.branchlength = 1e-7
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
    seqindextonodeindex = zeros(Int,length(sequences))
    seqindex = 1
    for (taxon,sequence) in zip(names, sequences)
        seqindextonodeindex[seqindex] = seqnametonodeindex[taxon]
        nodelist[seqnametonodeindex[taxon]].seqindex = seqindex
        seqindex += 1
    end

    markednodes = zeros(Int,length(nodelist))
    traitcounts = ones(Int, length(nodelist), 2)*-1
    stack = Int[1]
    while length(stack) > 0
        nodeindex = stack[end]
        currentnode = nodelist[nodeindex]
        if isleafnode(currentnode)
            trait = traits[seqnametoindex[currentnode.name]]
            for a=1:size(traitcounts,2)
                traitcounts[nodeindex,a] = 0
            end
            if trait > 0
                traitcounts[nodeindex, trait] = 1
            end
            pop!(stack)
        else
            cont = true
            for child in currentnode.children
                if traitcounts[child.nodeindex,1] == -1
                    push!(stack, child.nodeindex)
                    cont = false
                end
            end
            if cont
                for a=1:size(traitcounts,2)
                    traitcounts[nodeindex,a] = 0
                end
                for child in currentnode.children
                    traitcounts[nodeindex,:] += traitcounts[child.nodeindex,:]
                end
                pop!(stack)
            end
        end
    end
    for i=1:length(nodelist)
        nodeindex = nodelist[i].nodeindex
        parentnode = nodelist[i].parent
        parentindex = 0
        if !isnull(parentnode)
            parentindex = get(parentnode).nodeindex
        end
        if !(parentindex == 0 || maximum(traitcounts[parentindex,:]) == sum(traitcounts[parentindex,:])) && maximum(traitcounts[nodeindex,:]) == sum(traitcounts[nodeindex,:])
            markednodes[nodeindex] = argmax(traitcounts[nodeindex,:])
        end
    end

    numsamples = 100
    numcats = 20
    aadata,subcolumnrefs = getaadata(nodelist,seqindextonodeindex,sequences)
    initialmu = 0.4759392212135193
    initialshape = 4.851123709188913
    minf, minx = optimizealignmentlikelihood(nodelist, numcols, aadata, subcolumnrefs,initialmu,initialshape,numcats, 5)
    println("Finished",minx)
    maxll, sitelikelihoodsconditionals = aa_alignmentloglikelihood(nodelist, numcols, aadata, subcolumnrefs, minx[1], minx[2], numcats)
    rates = discretizegamma(minx[2], 1.0/minx[2], numcats)
    fout = open(outfile, "w")
    #,\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\"
    println(fout, "\"Site\",\"Target AA\",\"Chi2: X->notX vs X->X\",\"Fisher's: X->notX vs X->X\",\"Chi2: notX->notX vs notX->X\",\"Fisher's: notX->notX vs notX->X\",\"Chi2: notX vs X\",\"Fisher's: notX vs X\"")
    for col=1:numcols
        println(col,"\t", fastafile)
        samples = zeros(Int, numsamples, length(nodelist))
        for s=1:numsamples
            index = CommonUtils.sample(rng, sitelikelihoodsconditionals[col,:])
            logm = zeros(Float64, length(nodelist))
            aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)
            Q = minx[1]*LGmatrix*rates[index]
            sample = backwardssampling(rng, nodelist, aalikelihoods, Q, LGfreqs, 1)
            for z=1:length(nodelist)
                samples[s,z] = sample[1,z]
            end
        end
        for targetaa=1:20
            chi1,chi2,chi3,fisher1,fisher2,fisher3 = chisquaredtestmedians(nodelist, markednodes, traits, samples, targetaa, true)
            println(fout,col,",",aminoacids[targetaa],",",chi1,",",fisher1,",",chi2,",",fisher2,",",chi3,",",fisher3)
        end
    end
    close(fout)

    return outfile
end

function getchar(path::Array{Int,1}, time::Array{Float64}, t::Float64)
    c = path[1]
    for (p1, t1) in zip(path,time)
        if (t1 > t)
            break
        end
        c = p1
    end
    return c
end

function HAscore(seq::Array{Int,1})
    score = 0.0
    total = 0.0
    if aminoacids[seq[1]] == 'T'
        score += 21.2103485709695
    end
    total += 21.2103485709695

    if aminoacids[seq[2]] == 'N'
        score += 5.619267763
    end
    total += 5.619267763

    if aminoacids[seq[3]] == 'I'
        score += 9.235710963
    end
    total += 9.235710963

    if aminoacids[seq[4]] == 'I'
        score += 6.755671367
    end
    total += 6.755671367

    if aminoacids[seq[5]] == 'N'
        score += 2.581189838
    end
    total += 2.581189838
    return score/total
end

function getjsontree(nodelist::Array{TreeNode,1}, traits::Array{Int,1})
    json = Dict{Any,Any}()
    json["numnodes"] = length(nodelist)

    for node in nodelist
        json[node.nodeindex] = get(json, node.nodeindex, Dict{AbstractString,Any}())
        json[node.nodeindex]["nodeindex"] = node.nodeindex
        json[node.nodeindex]["branchlength"] = node.branchlength
        json[node.nodeindex]["nodedepth"] = getnodedepth(node)
        json[node.nodeindex]["rootnodedistance"] = rootnodedistance(node)
        json[node.nodeindex]["leafcount"] = countleafnodes(node)
        json[node.nodeindex]["name"] = node.name
        if isleafnode(node)
            if traits[node.seqindex] == 1
                json[node.nodeindex]["nodecolor"] = "#0000FF"
            elseif traits[node.seqindex] == 2
                json[node.nodeindex]["nodecolor"] = "#FF0000"
            end
        end
        if isnull(node.parent)
            json[node.nodeindex]["parentindex"] = -1
        else
            json[node.nodeindex]["parentindex"] = get(node.parent).nodeindex
        end
        json[node.nodeindex]["children"] = Int[child.nodeindex for child in node.children]
        childcount = 1
        for child in node.children
            json[child.nodeindex] = get(json, child.nodeindex, Dict{AbstractString,Any}())
            json[child.nodeindex]["childindex"] = childcount
            childcount += 1
        end
    end


    json[1]["lower"] = 0.0
    json[1]["upper"] = 1.0
    for node in nodelist
        lower = json[node.nodeindex]["lower"]
        upper = json[node.nodeindex]["upper"]
        json[node.nodeindex]["mid"] = lower +((upper-lower)/2.0)
        delta = upper-lower
        currenty = lower
        println("node: ", node.nodeindex)
        for child in node.children
            json[child.nodeindex]["lower"] = currenty
            json[child.nodeindex]["upper"] = currenty+(json[child.nodeindex]["leafcount"])/(json[node.nodeindex]["leafcount"])*delta
            currenty = json[child.nodeindex]["upper"]
            println(child.nodeindex,"\t",child.name,"\t",json[child.nodeindex]["leafcount"],"\t",json[node.nodeindex]["leafcount"],"\t",json[child.nodeindex]["lower"],"\t",json[child.nodeindex]["upper"])
        end
    end
    return JSON.json(json)
end

#main()
#mcmc()
#bhattacharya()
