using Pkg

println("Installing Julia packages...")
Pkg.add("Nullables")
Pkg.add("FastaIO")
Pkg.add("NLopt")
Pkg.add("Formatting")
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("ArgParse")
Pkg.add("JSON")
Pkg.add("HypothesisTests")
Pkg.add("Random")
Pkg.add("LinearAlgebra")
Pkg.add("SHA")
Pkg.add("Printf")
Pkg.add("CSV")
#Pkg.add("ROC")
Pkg.add("Plots")
Pkg.add("NPZ")
println("Running...")

module ParallelEvolution

push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using BranchPaths
using Binaries
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
using CSV

aminoacids = "ACDEFGHIKLMNPQRSTVWY"

function chi2pval(D, n)
        if D <= 0.0
                return 1.0
        end
        return ccdf(Chisq(n), D)
end

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
        "--collapsetraits"
        help = ""
        arg_type = String
        "--output"
        help = ""
        arg_type = String
        #default = "outfile.txt"
        "--mlmodel"
        help = ""
        action = :store_true
        "--mlindelmodel"
        help = ""
        action = :store_true
        "--score"
        help = ""
        action = :store_true
        "--coevolution"
        help = ""
        action = :store_true
        "--simulation"
        help = ""
        action = :store_true
        "--processsimulations"
        help = ""
        action = :store_true
    end

    return parse_args(s)
end

function getIndelTmatrix(sigma::Float64, delfreq::Float64, tau::Float64, hpfreq::Float64, lambda::Float64)
    indelfreqs = zeros(Float64,2)
    indelfreqs[1] = delfreq
    indelfreqs[2] = 1.0 - delfreq
    traitfreqs = zeros(Float64,2)
    traitfreqs[1] = hpfreq
    traitfreqs[2] = 1.0 - hpfreq

    DEL = 1
    INS = 2
    HP = 1
    LP = 2
    freqs = zeros(4)
    for b=1:2
        for n=1:2
            j = (b-1)*2 + n
            freqs[j] = indelfreqs[b]*traitfreqs[n]
            if n == HP
                if b == DEL
                    freqs[j] *= lambda
                else
                    freqs[j] /= lambda
                end
            end
        end
    end
    freqs /= sum(freqs)

    Q = zeros(4,4)
    for a=1:2
        for b=1:2
            for m=1:2
                for n=1:2
                    i = (a-1)*2 + m
                    j = (b-1)*2 + n
                    
                    if (a != b && m == n) || (a == b && m != n)
                        if a != b
                            rate = indelfreqs[b]*sigma
                        end
                        if m != n
                            rate = traitfreqs[n]*tau
                        end

                        if m == n && n == HP
                            if b == DEL
                                rate *= lambda
                            elseif b == INS
                                rate /= lambda
                            end
                        elseif a == b && n == HP
                            if b == DEL
                                rate *= lambda
                            elseif b == INS
                                rate /= lambda
                            end
                        end
                        Q[i,j] = rate
                        Q[i,i] -= rate
                    end
                end
            end
        end
    end

    return Q, freqs
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

function getindelcolumn(sequences::Array{AbstractString,1}, col::Int)
    ret = zeros(Int, length(sequences))
    s = 1
    for seq in sequences
        if seq[col] == '-'
            ret[s] = 0
        else 
            ret[s] = 1
        end
        s += 1
    end
    return ret
end



function backstack(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, logm::Array{Float64,1}, initialfreqs::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, alphabet::Int)
    backlikelihoods = copy(likelihoods)
    backlikelihoods[1,:] = backlikelihoods[1,:].*initialfreqs
    for node in nodelist[2:end]
        parentindex = get(node.parent).nodeindex
        backlikelihoods[node.nodeindex,:] = backlikelihoods[parentindex,:].*(transprobs[node.nodeindex]*likelihoods[node.nodeindex,:])
        backlikelihoods[node.nodeindex,:] /= sum(backlikelihoods[node.nodeindex,:])
    end
    return backlikelihoods
end

function marginalfreqs(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, logm::Array{Float64,1}, initialfreqs::Array{Float64,1}, transprobs::Array{Array{Float64,2},1}, alphabet::Int)
    forwardliks = felsensteinstack(nodelist, likelihoods, logm, transprobs, alphabet)
    backliks = backstack(nodelist, forwardliks, logm, initialfreqs, transprobs, alphabet)    
    marginalfreqs = copy(forwardliks)
    marginalfreqs[1,:] /= sum(marginalfreqs[1,:])
    for node in nodelist[2:end]
        marginalfreqs[node.nodeindex,:] = forwardliks[node.nodeindex,:].*backliks[node.nodeindex,:]
        marginalfreqs[node.nodeindex,:] /= sum(marginalfreqs[node.nodeindex,:])
    end
    return marginalfreqs
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
    #println("freqs",freqs)
    #println("Q",Q)
    #println(transprobs[2])
    #println("blength", nodelist[2].branchlength)

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
                    v = transprobs[child.nodeindex][parentc,:].*likelihoods[child.nodeindex,:]
                    z = CommonUtils.sample(rng, v)
                    #=
                    if length(v) == 40
                        println(child.nodeindex,"\t",v,"\t",v[parentc],"\t",parentc,"\t",z)
                    end=#
                    samples[s, child.nodeindex] = z
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
        println(e)
        exit()
        return -1e20
    end
end


function optimizeaatraitmodel(datalikelihoods::Array{Float64,2},mles::Array{Float64,2}, mlparams::Dict{Tuple{Int,Int}, AATraitParameters}, nodelist::Array{TreeNode,1}, data::Array{Float64,3}, subcolumnrefs::Array{Int,2}, initialparams::AATraitParameters, targetaa::Int, col::Int, nloptmethod::Symbol=:LN_COBYLA, fixtraitparams::Bool=false)
    mode = 1
    initialarr = getarray(initialparams)
    #println("A",initialarr)
    global optiters
    optiters = 0
    maxoptiter = 1000
    opt = Opt(nloptmethod, 4)
    #localObjectiveFunction = ((param, grad) -> computesiteloglikelihoods(mles, mlparams, nodelist, data, subcolumnrefs, AATraitParameters(param), targetaa)[col])
    localObjectiveFunction = ((param, grad) -> computesiteloglikelihoodhelper(param, datalikelihoods, mles, mlparams, nodelist, data, subcolumnrefs, initialparams, targetaa, col, 100))

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
        xtol_rel!(opt,1e-4)
        maxeval!(opt, maxoptiter)
    else
        xtol_rel!(opt,1e-6)
        maxeval!(opt, 3*maxoptiter)
    end
    max_objective!(opt, localObjectiveFunction)
    #println("B",initialarr)
    initialarr[4] = max(0.011,initialarr[4])
    initialarr[4] = min(0.989,initialarr[4])
    (maxll,maxparams,ret) = optimize(opt, initialarr)
    println("Optimal: ", maxll,"\t",maxparams)
    return maxll,maxparams
end

function optimizeaamodel(aalikelihoods::Array{Float64,2},mles::Array{Float64,2}, mlparams::Dict{Tuple{Int,Int}, AATraitParameters}, nodelist::Array{TreeNode,1}, data::Array{Float64,3}, subcolumnrefs::Array{Int,2}, initialparams::AATraitParameters, targetaa::Int, col::Int, nloptmethod::Symbol=:LN_COBYLA, fixtraitparams::Bool=false)
    mode = 1
    initialarr = Float64[1.0]
    global optiters
    optiters = 0
    maxoptiter = 2000
    opt = Opt(nloptmethod, 1)

    localObjectiveFunction = ((param, grad) -> getaacolumnloglikelihood(nodelist, aalikelihoods, param[1]))

    lower = ones(Float64, 1)*1e-3
    upper = ones(Float64, 1)*1e3
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

function getindeltraitlikelihoods(nodelist::Array{TreeNode,1}, sequences::Array{AbstractString,1}, aacolumn::Array{Int,1}, traitcolumn::Array{Int,1}, seqindextonodeindex::Array{Int,1})
    alphabet = 4
    likelihoods = ones(length(nodelist),alphabet)*-Inf
    for seqindex=1:length(sequences)
        nodeindex = seqindextonodeindex[seqindex]
        indel = 1
        if aacolumn[seqindex] != 0
            indel = 2
        end
        for a=1:alphabet
            likelihoods[nodeindex,a] = 0.0
        end
        if traitcolumn[seqindex] == 0 && indel == 0
            for a=1:alphabet
                likelihoods[nodeindex,a] = 1.0
            end
        elseif traitcolumn[seqindex] == 0
            likelihoods[nodeindex,indel] = 1.0
            likelihoods[nodeindex,2+indel] = 1.0
        elseif  indel == 0
            for a=1:2
                likelihoods[nodeindex,(traitcolumn[seqindex]-1)*2+a] = 1.0
            end
        else
            likelihoods[nodeindex,(traitcolumn[seqindex]-1)*2 + indel] = 1.0
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

function secondderiv(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, params::AATraitParameters, targetaa::Int)
    maxparams = deepcopy(params)
    lambdas = Float64[]
    liks = Float64[]
    cumll = -Inf
    cumliks = Float64[]
    mid = 500
    v = exp(log(100.0)/mid)
    for i=1:1001
        index = i - mid - 1
        multiplier = v^index
        #push!(lambdas, maxparams.lambda*multiplier)
        push!(lambdas, v*multiplier)
        tempparams = deepcopy(maxparams)
        tempparams.lambda = lambdas[end]
        ll = -Inf
        try
            ll = getaatraitcolumnloglikelihood(nodelist, likelihoods, tempparams, targetaa)
            if i == 1
                ll -= log(lambdas[end]-0.0)
            else
                ll -= log(lambdas[end]-lambdas[end-1])
            end
        catch

        end
        push!(liks,ll)
        cumll = logsumexp(cumll, ll)
        push!(cumliks,cumll)        
      
    end
    cumliks = exp.(cumliks .- maximum(cumliks))
    lower = lambdas[1]
    upper = lambdas[end]
    for i=1:1000
        if cumliks[i] <= 0.025
            lower = lambdas[i]
        elseif cumliks[i] >= 0.975
            upper = lambdas[i]
            break
        end
    end
    #println(lambdas)
    println((lower,upper,1.0-cumliks[mid+1]))
    #=   
    maxparams = deepcopy(params)
    h = maxparams.lambda*0.01
    lambdas = Float64[maxparams.lambda-h, maxparams.lambda, maxparams.lambda+h]
    liks = Float64[]
    for lambda in lambdas
        tempparams = deepcopy(maxparams)
        tempparams.lambda = lambda
        ll = getaatraitcolumnloglikelihood(nodelist, likelihoods, tempparams, targetaa)
        push!(liks,ll)
    end
    deriv = (liks[3]-2.0*liks[2]+liks[1])/(h*h)
    interval = 1.96/sqrt(-deriv)
    println("lambdas: ", lambdas)
    println("liks: ", liks)
    println("sec: ", deriv,"\t", interval)
    return interval=#
    return (lower,upper,1.0-cumliks[mid+1])
end

function confintervalsindeltrait(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, params::Array{Float64,1})
    maxparams = deepcopy(params)
    lambdas = Float64[]
    liks = Float64[]
    cumll = -Inf
    cumliks = Float64[]
    mid = 500
    v = exp(log(100.0)/mid)
    for i=1:1001
        index = i - mid - 1
        multiplier = v^index
        #push!(lambdas, maxparams.lambda*multiplier)
        push!(lambdas, v*multiplier)
        tempparams = deepcopy(maxparams)
        tempparams[5] = lambdas[end]
        ll = -Inf
        try
            ll = getindeltraitcolumnloglikelihood(nodelist, likelihoods, tempparams)
            if i == 1
                ll -= log(lambdas[end]-0.0)
            else
                ll -= log(lambdas[end]-lambdas[end-1])
            end
        catch

        end
        push!(liks,ll)
        cumll = logsumexp(cumll, ll)
        push!(cumliks,cumll)        
      
    end
    cumliks = exp.(cumliks .- maximum(cumliks))
    lower = lambdas[1]
    upper = lambdas[end]
    for i=1:1000
        if cumliks[i] <= 0.025
            lower = lambdas[i]
        elseif cumliks[i] >= 0.975
            upper = lambdas[i]
            break
        end
    end
    #println(lambdas)
    #println((lower,upper,1.0-cumliks[mid+1]))
    #=   
    maxparams = deepcopy(params)
    h = maxparams.lambda*0.01
    lambdas = Float64[maxparams.lambda-h, maxparams.lambda, maxparams.lambda+h]
    liks = Float64[]
    for lambda in lambdas
        tempparams = deepcopy(maxparams)
        tempparams.lambda = lambda
        ll = getaatraitcolumnloglikelihood(nodelist, likelihoods, tempparams, targetaa)
        push!(liks,ll)
    end
    deriv = (liks[3]-2.0*liks[2]+liks[1])/(h*h)
    interval = 1.96/sqrt(-deriv)
    println("lambdas: ", lambdas)
    println("liks: ", liks)
    println("sec: ", deriv,"\t", interval)
    return interval=#
    return (lower,upper,1.0-cumliks[mid+1])
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

function getindeltraitcolumnloglikelihood(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, params::Array{Float64,1})
    alphabet = 4
    Q,freqs = getIndelTmatrix(params[1],params[2],params[3],params[4],params[5])
    logm = zeros(Float64, size(likelihoods,1))
    felsensteinstack(nodelist, likelihoods,logm, gettransprobmatrices(nodelist, Q), alphabet)
    totalloglikelihood = logm[1]+log.(dot(likelihoods[1,:], freqs))
    return totalloglikelihood
end

function optimizeindeltraitmodel(datalikelihoods::Array{Float64,2}, nodelist::Array{TreeNode,1}, nloptmethod::Symbol=:LN_COBYLA; fixlambda::Bool=false, initialparams::Array{Float64,1}=Float64[1.0, 0.2, 1.0, 0.2, 1.0])
    global optiters
    optiters = 0
    maxoptiter = 1000
    opt = Opt(nloptmethod, 5)
    localObjectiveFunction = ((param, grad) -> getindeltraitcolumnloglikelihood(nodelist, datalikelihoods, param))

    lower = ones(Float64, 5)*1e-3
    upper = ones(Float64, 5)*1e3
    lower[2] = 0.0
    lower[4] = 0.0
    upper[2] = 1.0
    upper[4] = 1.0
    if fixlambda
        lower[5] = 1.0
        upper[5] = 1.0
    end
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)

    maxeval!(opt, 2000)
    max_objective!(opt, localObjectiveFunction)
    (maxll,maxparams,ret) = optimize(opt, initialparams)
    return maxll,maxparams
end

function getaacolumnloglikelihood(nodelist::Array{TreeNode,1}, likelihoods::Array{Float64,2}, mu::Float64)
    logm = zeros(Float64, size(likelihoods,1))
    felsensteinstack(nodelist, likelihoods,logm, gettransprobmatrices(nodelist, mu*LGmatrix), 20)
    totalloglikelihood = logm[1]+log.(dot(likelihoods[1,:], LGfreqs))
    #println("getaacolumnloglikelihood ",mu,"\t",totalloglikelihood)
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

function mlmodel()
    parsed_args = parse_commandline()
    fastafile = parsed_args["alignment"]
    treefile = parsed_args["tree"]
    outfilename = parsed_args["output"]
    annotations = parsed_args["annotations"]
    return mlmodel(MersenneTwister(1),fastafile,treefile,outfilename,annotations)
end

function mlmodel(rng,fastafile,treefile,outfilename,annotations,maxcols=10000000000,prioritycolumns=Int[])
    parsed_args = parse_commandline()
    if outfilename == nothing
        outfilename = string(fastafile,".results.csv")
    end
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
        end
    end
    #println(sequences)
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
    annotationnames = reverse(split(annotations,","))

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
    #println(traits)

    nodelist = loadtree(treefile)
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


    #performchi2association(string(fastafile,".chi2.csv"), sequences, traits, annotationnames)

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
    for c=1:min(numcols,maxcols)
        if !in(c,columns)
            push!(columns,c)
        end
    end

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
    confintlower = zeros(numcols,21)
    confintupper = zeros(numcols,21)
    posterior = zeros(numcols,21)
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

    #println("a", tll1)
    #println("b", tll2)
    p2 = AATraitParameters(Float64[31.47925915,1.0,6.040222915,0.24263616])
    p1 = AATraitParameters(Float64[31.47925915,1.0,7.456053451, 0.394240064])
    cll1 = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, p2, 0)
    cll2 = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, p1, 0)
    aall = getaacolumnloglikelihood(nodelist,aalikelihoods, 31.47925915)
    #println("e",aall)
    #println("c",cll1,"\t", tll1+aall, "\t",tll1+aall-cll1)
    #println("d",cll2,"\t", tll2+aall, "\t",tll2+aall-cll2)
    #println(tll1-tll2)
    #println(cll1-cll2)

    #initial = AATraitParameters(initialparams)
    initial = deepcopy(p1)
    #println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_DIRECT)
    #println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_CRS2_LM)
    #println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :LN_COBYLA)
    #println(getarray(initial))
    initialparams = getarray(initial)

    
    marginalsfile = string(outfilename,".marginals")
    outmarginals = open(marginalsfile,"w")     
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
                    if targetaa == 0
                        #println("COL $(col)")
                        temp = getaacolumn(sequences, col)
                        aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)
                        maxll,maxparams = optimizeaamodel(aalikelihoods,mles, mlparams, nodelist, data, subcolumnrefs, initial, targetaa, col, nloptmethod, targetaa == 0)
                        initial.mu = maxparams[1]
                        mlparams[key] = deepcopy(initial)
                        mles[col,21] = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, initial, 0)
                        #println("HERE")
                        #println(mlparams[key])
                        #println(mles[col,21])
                    end

                    maxll,maxparams = optimizeaatraitmodel(datalikelihoods,mles, mlparams, nodelist, data, subcolumnrefs, initial, targetaa, col, nloptmethod, targetaa == 0)
                    mlparams[key] = AATraitParameters(maxparams)
                    Q, freqs, marginaltraitfreqs = getAATraitmatrix(AATraitParameters(maxparams), LGmatrix, LGfreqs, targetaa)
                    #println(marginaltraitfreqs)
                    if targetaa == 0
                        mles[col,21] = maxll
                        optimisationstatus[col,21] = 1
                        initialparams[2:4] = maxparams[2:4]
                        samples = backwardssampling(rng,nodelist, datalikelihoods, Q, freqs, 100)

                        #ZZZZZZZ
                        aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)                    
                        logm = zeros(Float64, size(aalikelihoods,1))
                        mu = initialparams[1]
                        marginals = marginalfreqs(nodelist, aalikelihoods, logm, LGfreqs, gettransprobmatrices(nodelist, mu*LGmatrix), 20)
                        targetmarginals = zeros(Float64,2,20)                        
                        for node in nodelist
                            if markednodes[node.nodeindex] > 0
                                maxfreq,maxaa = findmax(marginals[node.nodeindex,:])
                                targetmarginals[markednodes[node.nodeindex], :] += marginals[node.nodeindex,:]
                                println(outmarginals, col,"\t",node.nodeindex,"\t", markednodes[node.nodeindex],"\t", aminoacids[maxaa],"\t",marginals[node.nodeindex,:])
                                flush(outmarginals)
                            end
                        end
                        targetmarginals[1,:] /= sum(targetmarginals[1,:])
                        targetmarginals[2,:] /= sum(targetmarginals[2,:])
                        println(outmarginals,"HP ", targetmarginals[1,:])
                        println(outmarginals,"LP ", targetmarginals[2,:])
                        print(outmarginals,"HP: ")
                        for t=1:20
                            if targetmarginals[1,t] > 0.001
                                print(outmarginals,aminoacids[t], " (",@sprintf("%0.4f", targetmarginals[1,t]),");  ")
                            end
                        end
                        println(outmarginals,"")
                        print(outmarginals,"LP: ")
                        for t=1:20
                            if targetmarginals[2,t] > 0.001
                                print(outmarginals,aminoacids[t], " (",@sprintf("%0.4f", targetmarginals[2,t]),");  ")
                            end
                        end
                        println(outmarginals,"")

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
                        #println("chi2=",@sprintf("%0.2f", D),"\tpval=",ccdf(Chisq(1), D))

                        #println("AAA: ",mlparams[(col,targetaa)])
                        confintlower[col,targetaa], confintupper[col,targetaa], posterior[col,targetaa] = secondderiv(nodelist,datalikelihoods,mlparams[(col,targetaa)],targetaa)
                    end
                    if targetaa == 0
                        initialparams = copy(maxparams)
                    end
                end
            end

            #println(abspath(outfilename))
            #exit()
            outfile = open(outfilename,"w")       
            #println(outfile,"\"Site\",\"Target AA\",\"Count at leaf\",\"mu\",\"tau\",\"p\",\"lambda\",\"Optimisation\",\"MLL\",\"AIC\",\"chi2\",\"pvalue\",\"X->notX vs X->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"X->notX vs X->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\",\"Fisher's pvalue\",\"lower\",\"upper\"")
            println(outfile,"\"Site\",\"Target AA\",\"Count at leaf\",\"mu\",\"tau\",\"p\",\"lambda\",\"lower 2.5%\",\"upper 2.5%\",\"posterior\",\"Optimisation\",\"MLL\",\"AIC\",\"chi2\",\"pvalue\",\"X->notX vs X->X\",\"pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"pvalue\",\"lower\",\"upper\",\"notX vs X\",\"pvalue\",\"lower\",\"upper\",\"X->notX vs X->X\",\"pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"pvalue\",\"lower\",\"upper\",\"notX vs X\",\"pvalue\",\"lower\",\"upper\"")
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
                               # println(outfile,c,",", "independent,","-,", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel, ",",mles[c,21], ",",aic, ",", "0.0", ",", "1.0", ",0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0")
                                println(outfile,c,",", "independent,","-,", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda,",-",",-",",0.5", ",",optstatuslabel, ",",mles[c,21], ",",aic, ",", "0.0", ",", "1.0", ",0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0")
                            else
                                optstatuslabel = "Partial"
                                if optimisationstatus[c,t] == 1
                                    optstatuslabel = "Complete"
                                end
                                D = 2.0*(mles[c,t]-mles[c,21])
                                pval = chi2pval(D,1)
                                if params.lambda < 1.0
                                    D = -D
                                end
                                aic = 2.0*4.0 - 2.0*mles[c,t]
                                #println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", fishertable_hponly[col,t,1], ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", fishertable_hponly[col,t,2], ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", fishertable_hponly[col,t,3], ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", fishertable[col,t,1], ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", fishertable[col,t,2], ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", fishertable[col,t,3], ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                
                               # println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,1])), ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,2])), ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,3])), ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians[c,t,1])), ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians[c,t,2])), ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians[c,t,3])), ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                if confintlower[c,t] == confintupper[c,t] == posterior[c,t] == 0.0
                                    println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda,",","-",",","-",",","0.5", ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,1])), ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,2])), ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,3])), ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians[c,t,1])), ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians[c,t,2])), ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians[c,t,3])), ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                else
                                    println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda,",",confintlower[c,t],",",confintupper[c,t],",",posterior[c,t], ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,1])), ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,2])), ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,3])), ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians[c,t,1])), ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians[c,t,2])), ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians[c,t,3])), ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                end
                            end
                        end
                    end
                end
            end
            close(outfile)
        end
    end
end

function mlindelmodel()
    parsed_args = parse_commandline()
    fastafile = parsed_args["alignment"]
    treefile = parsed_args["tree"]
    outfilename = parsed_args["output"]
    annotations = parsed_args["annotations"]
    return mlindelmodel(MersenneTwister(1),fastafile,treefile,outfilename,annotations)
end

function mlindelmodel(rng,fastafile,treefile,outfilename,annotations,maxcols=10000000000,prioritycolumns=Int[])
    parsed_args = parse_commandline()
    if outfilename == nothing
        outfilename = string(fastafile,".results.update.csv")
    end
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
        end
    end
    #println(sequences)
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
    annotationnames = reverse(split(annotations,","))

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
    #println(traits)

    nodelist = loadtree(treefile)
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


    #performchi2association(string(fastafile,".chi2.csv"), sequences, traits, annotationnames)

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
    for c=1:min(numcols,maxcols)
        if !in(c,columns)
            push!(columns,c)
        end
    end

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
    confintlower = zeros(numcols,21)
    confintupper = zeros(numcols,21)
    posterior = zeros(numcols,21)
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

    startindex = 1
    endindex = length(sequences[1])
    startindex = 145
    endindex = 145
    for col=startindex:endindex
        initialparams = Float64[1.0,0.2,1.0,0.2,1.0]
        if in(0, getindelcolumn(sequences, col))      
            datalikelihoods = getindeltraitlikelihoods(nodelist, sequences, getindelcolumn(sequences, col), traits, seqindextonodeindex)
            likelihood = getindeltraitcolumnloglikelihood(nodelist, datalikelihoods, initialparams)
            maxll1, initialparams = optimizeindeltraitmodel(datalikelihoods, nodelist,fixlambda=true,initialparams=initialparams)
            maxll2, maxparams = optimizeindeltraitmodel(datalikelihoods, nodelist,fixlambda=false,initialparams=initialparams)
            lower,upper,posterior = confintervalsindeltrait(nodelist,datalikelihoods,maxparams)
            println(col,"\t", chi2pval(2.0*(maxll2-maxll1),1), "\t",maxll1,"\t",maxll2,"\t",maxparams[5],"\t",lower,"\t",upper,"\t",posterior)
            Q,freqs = getIndelTmatrix(maxparams[1], maxparams[2], maxparams[3], maxparams[4], maxparams[5])        
        end
    end
end

function score()
    parsed_args = parse_commandline()
    fastafile = parsed_args["alignment"]
    treefile = parsed_args["tree"]
    outfilename = parsed_args["output"]
    annotations = parsed_args["annotations"]
    return score(MersenneTwister(1),fastafile,treefile,outfilename,annotations)
end

function score(rng,fastafile,treefile,outfilename,annotations,maxcols=10000000000,prioritycolumns=Int[])
    parsed_args = parse_commandline()

    if outfilename == nothing
        outfilename = string(fastafile,".results.update.csv")
    end
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
        end
    end
    println(sequences)
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
    annotationnames = reverse(split(annotations,","))

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

    nodelist = loadtree(treefile)
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


    #performchi2association(string(fastafile,".chi2.csv"), sequences, traits, annotationnames)

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
    for c=1:min(numcols,maxcols)
        if !in(c,columns)
            push!(columns,c)
        end
    end

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
    confintlower = zeros(numcols,21)
    confintupper = zeros(numcols,21)
    posterior = zeros(numcols,21)
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

    logm = zeros(Float64, size(aalikelihoods,1))
    mu = initialparams[1]
    freqs = marginalfreqs(nodelist, aalikelihoods, logm, LGfreqs, gettransprobmatrices(nodelist, mu*LGmatrix), 20)
    for node in nodelist
        println(node.nodeindex,"\t",freqs[node.nodeindex,:])
    end
    
    #felsensteinstack(nodelist, aalikelihoods,logm, gettransprobmatrices(nodelist, mu*LGmatrix), 20)
    #totalloglikelihood = logm[1]+log.(dot(aalikelihoods[1,:], LGfreqs))
    #println(aalikelihoods[1,:]/sum(aalikelihoods[1,:]))

    #=
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
                    if targetaa == 0
                        #println("COL $(col)")
                        temp = getaacolumn(sequences, col)
                        aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)
                        maxll,maxparams = optimizeaamodel(aalikelihoods,mles, mlparams, nodelist, data, subcolumnrefs, initial, targetaa, col, nloptmethod, targetaa == 0)
                        initial.mu = maxparams[1]
                        mlparams[key] = deepcopy(initial)
                        mles[col,21] = getaatraitcolumnloglikelihood(nodelist, datalikelihoods, initial, 0)
                        #println("HERE")
                        #println(mlparams[key])
                        #println(mles[col,21])
                    end

                    maxll,maxparams = optimizeaatraitmodel(datalikelihoods,mles, mlparams, nodelist, data, subcolumnrefs, initial, targetaa, col, nloptmethod, targetaa == 0)
                    mlparams[key] = AATraitParameters(maxparams)
                    Q, freqs, marginaltraitfreqs = getAATraitmatrix(AATraitParameters(maxparams), LGmatrix, LGfreqs, targetaa)
                    #println(marginaltraitfreqs)
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
                        #println("chi2=",@sprintf("%0.2f", D),"\tpval=",ccdf(Chisq(1), D))

                        println("AAA: ",mlparams[(col,targetaa)])
                        confintlower[col,targetaa], confintupper[col,targetaa], posterior[col,targetaa] = secondderiv(nodelist,datalikelihoods,mlparams[(col,targetaa)],targetaa)
                    end
                    if targetaa == 0
                        initialparams = copy(maxparams)
                    end
                end
            end

            #println(abspath(outfilename))
            #exit()
            outfile = open(outfilename,"w")
            #println(outfile,"\"Site\",\"Target AA\",\"Count at leaf\",\"mu\",\"tau\",\"p\",\"lambda\",\"Optimisation\",\"MLL\",\"AIC\",\"chi2\",\"pvalue\",\"X->notX vs X->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"X->notX vs X->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"Fisher's pvalue\",\"lower\",\"upper\",\"notX vs X\",\"Fisher's pvalue\",\"lower\",\"upper\"")
            println(outfile,"\"Site\",\"Target AA\",\"Count at leaf\",\"mu\",\"tau\",\"p\",\"lambda\",\"lower 2.5%\",\"upper 2.5%\",\"posterior\",\"Optimisation\",\"MLL\",\"AIC\",\"chi2\",\"pvalue\",\"X->notX vs X->X\",\"pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"pvalue\",\"lower\",\"upper\",\"notX vs X\",\"pvalue\",\"lower\",\"upper\",\"X->notX vs X->X\",\"pvalue\",\"lower\",\"upper\",\"notX->notX vs notX->X\",\"pvalue\",\"lower\",\"upper\",\"notX vs X\",\"pvalue\",\"lower\",\"upper\"")
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
                               # println(outfile,c,",", "independent,","-,", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel, ",",mles[c,21], ",",aic, ",", "0.0", ",", "1.0", ",0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0")
                                println(outfile,c,",", "independent,","-,", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda,",-",",-",",0.5", ",",optstatuslabel, ",",mles[c,21], ",",aic, ",", "0.0", ",", "1.0", ",0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0")
                            else
                                optstatuslabel = "Partial"
                                if optimisationstatus[c,t] == 1
                                    optstatuslabel = "Complete"
                                end
                                D = 2.0*(mles[c,t]-mles[c,21])
                                pval = chi2pval(D,1)
                                if params.lambda < 1.0
                                    D = -D
                                end
                                aic = 2.0*4.0 - 2.0*mles[c,t]
                                #println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", fishertable_hponly[col,t,1], ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", fishertable_hponly[col,t,2], ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", fishertable_hponly[col,t,3], ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", fishertable[col,t,1], ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", fishertable[col,t,2], ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", fishertable[col,t,3], ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                
                               # println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda, ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,1])), ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,2])), ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,3])), ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians[c,t,1])), ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians[c,t,2])), ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians[c,t,3])), ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                if confintlower[c,t] == confintupper[c,t] == posterior[c,t] == 0.0
                                    println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda,",","-",",","-",",","0.5", ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,1])), ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,2])), ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,3])), ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians[c,t,1])), ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians[c,t,2])), ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians[c,t,3])), ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                else
                                    println(outfile,c,",", aminoacids[t],",",aminoacidatleaf[c,t],",", params.mu, ",",  params.tau, ",",  params.p, ",", params.lambda,",",confintlower[c,t],",",confintupper[c,t],",",posterior[c,t], ",",optstatuslabel,",",mles[c,t],",",aic, ",", D, ",", pval, ",", chi2medians_hponly[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,1])), ",", chi2lower_hponly[c,t,1], ",", chi2upper_hponly[c,t,1], ",", chi2medians_hponly[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,2])), ",", chi2lower_hponly[c,t,2], ",", chi2upper_hponly[c,t,2], ",", chi2medians_hponly[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians_hponly[c,t,3])), ",", chi2lower_hponly[c,t,3], ",", chi2upper_hponly[c,t,3], ",", chi2medians[c,t,1], ",", ccdf(Chisq(1), abs(chi2medians[c,t,1])), ",", chi2lower[c,t,1], ",", chi2upper[c,t,1], ",", chi2medians[c,t,2], ",", ccdf(Chisq(1), abs(chi2medians[c,t,2])), ",", chi2lower[c,t,2], ",", chi2upper[c,t,2], ",", chi2medians[c,t,3], ",", ccdf(Chisq(1), abs(chi2medians[c,t,3])), ",", chi2lower[c,t,3], ",", chi2upper[c,t,3])
                                end
                            end
                        end
                    end
                end
            end
            close(outfile)
        end
    end=#
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
            push!(sequences, replace(seq, "\r" => ""))
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
    annotationnames = reverse(split(parsed_args["annotations"],","))

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

    nodelist = loadtree(treefile)
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

function getLGmatrixwithdeletion(deletionrate::Float64)
    Q = zeros(Float64,21,21)
    for aa1=1:20
        for aa2=1:20
            if aa1 != aa2
                Q[aa1,aa2] = LGexchangeability[aa1,aa2]*LGfreqs[aa2]
            end
        end
    end
    for aa=1:20
        Q[21,aa] = deletionrate*LGfreqs[aa]
        Q[aa,21] = Q[21,aa]
    end
    for aa1=1:21
        Q[aa1,aa1] = 0.0
        for aa2=1:21
            if aa1 != aa2
                Q[aa1,aa1] -= Q[aa1,aa2]
            end
        end
    end
    return Q
end

function getLGmatrixwithdeletion(deletionfreq::Float64, deletionrate::Float64)
    S = zeros(Float64,21,21)
    freqs = zeros(Float64,21)
    for aa1=1:20
        freqs[aa1] = LGfreqs[aa1]
        for aa2=1:20
            S[aa1,aa2] = LGexchangeability[aa1,aa2]
        end
    end 
    for aa=1:20
        S[21,aa] = deletionrate
        S[aa,21] = deletionrate
    end
    freqs *= (1.0-deletionfreq)
    freqs[21] = deletionfreq


    Q = zeros(Float64,21,21)
    for aa1=1:21
        for aa2=1:21
            if aa1 != aa2
                Q[aa1,aa2] = S[aa1,aa2]*freqs[aa2]
            end
        end
    end

    for aa1=1:21
        Q[aa1,aa1] = 0.0
        for aa2=1:21
            if aa1 != aa2
                Q[aa1,aa1] -= Q[aa1,aa2]
            end
        end
    end

    println(S,"\t",freqs,"\t",Q)
    return Q,freqs
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
    collapsetraits = Int[]
    if parsed_args["collapsetraits"] != nothing
        collapsetraits = Int[parse(Int,s) for s in split(parsed_args["collapsetraits"],",")]
    end
    return bhattacharya(parsed_args["alignment"], parsed_args["tree"], split(parsed_args["annotations"],","), collapsetraits, parsed_args["output"])
end

function bhattacharya(fastafile::AbstractString, treefile, annotationnames, collapsetraits::Array{Int,1}, outfile::AbstractString)
    println("COLLAPSED TRAITS ", collapsetraits)
    if treefile == nothing
        treestring, treefile = Binaries.fasttreeaa(fastafile)
    end

    rng = MersenneTwister(2184104820249809)
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
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

    nodelist = loadtree(treefile)

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
println(seqnametonodeindex)
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

#collapsetraits = Int[1,2]
nodetraits = zeros(Int, length(nodelist))
for i=1:length(traits)
    if !(traits[i] in collapsetraits)
        nodetraits[seqindextonodeindex[i]] = traits[i]
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
        traitvalue =  argmax(traitcounts[nodeindex,:])
        markednodes[nodeindex] = traitvalue
        if traitvalue in collapsetraits
            nodetraits[nodeindex] = traitvalue
        end
    end
end
#println("markednodes ", markednodes)
#println("nodetraits  ", nodetraits)

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
        chi1,chi2,chi3,fisher1,fisher2,fisher3 = chisquaredtestmedians(nodelist, nodetraits, samples, targetaa)
        #chi1,chi2,chi3,fisher1,fisher2,fisher3 = chisquaredtestmedians(nodelist, markednodes, traits, samples, targetaa, true)
        println(fout,col,",",aminoacids[targetaa],",",chi1,",",fisher1,",",chi2,",",fisher2,",",chi3,",",fisher3)
    end
end
close(fout)

return outfile
end

function coevolution()
    parsed_args = parse_commandline()
    collapsetraits = Int[]
    if parsed_args["collapsetraits"] != nothing
        collapsetraits = Int[parse(Int,s) for s in split(parsed_args["collapsetraits"],",")]
    end
    return coevolution(parsed_args["alignment"], parsed_args["tree"], split(parsed_args["annotations"],","), collapsetraits, parsed_args["output"])
end

function coevolution(fastafile::AbstractString, treefile, annotationnames, collapsetraits::Array{Int,1}, outfile::AbstractString)
    println("COLLAPSED TRAITS ", collapsetraits)
    if treefile == nothing
        treestring, treefile = Binaries.fasttreeaa(fastafile)
    end

    rng = MersenneTwister(2184104820249809)
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
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

    nodelist = loadtree(treefile)
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
    println(seqnametonodeindex)
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

    #collapsetraits = Int[1,2]
    nodetraits = zeros(Int, length(nodelist))
    for i=1:length(traits)
        if !(traits[i] in collapsetraits)
            nodetraits[seqindextonodeindex[i]] = traits[i]
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
            traitvalue =  argmax(traitcounts[nodeindex,:])
            markednodes[nodeindex] = traitvalue
            if traitvalue in collapsetraits
                nodetraits[nodeindex] = traitvalue
            end
        end
    end
    #println("markednodes ", markednodes)
    #println("nodetraits  ", nodetraits)

    initial = AATraitParameters(Float64[31.47925915,1.0,7.456053451, 0.394240064])
    traitlikelihoods = gettraitlikelihoods(nodelist, traits, seqindextonodeindex)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_DIRECT)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_CRS2_LM)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :LN_COBYLA)
    println(getarray(initial))

    numsamples = 250
    numcats = 20
    aadata,subcolumnrefs = getaadata(nodelist,seqindextonodeindex,sequences)
    initialmu = 0.45005877688975415
    initialshape = 3.8950017807139856
    minf, minx = optimizealignmentlikelihood(nodelist, numcols, aadata, subcolumnrefs,initialmu,initialshape,numcats, 1)
    println("Finished",minx)
    maxll, sitelikelihoodsconditionals = aa_alignmentloglikelihood(nodelist, numcols, aadata, subcolumnrefs, minx[1], minx[2], numcats)
    rates = discretizegamma(minx[2], 1.0/minx[2], numcats)
    #=
    samples = zeros(Int, col, numsamples, length(nodelist))

    sampledrateindices = zeros(Int,numcols)
    for z=1:numcols
        sampledrateindices[col] = CommonUtils.sample(rng, sitelikelihoodsconditionals[col,:])
    end=#


    selectedcols = []
    #=
    push!(selectedcols, (143, "T"))
    push!(selectedcols, (274, "I"))
    push!(selectedcols, (438, "I"))
    push!(selectedcols, (384, "N"))
    =#
    push!(selectedcols, (136, "N"))
    push!(selectedcols, (143, "T"))
    #push!(selectedcols, (402, "K"))
    #push!(selectedcols, (335, "T"))
    push!(selectedcols, (274, "I"))
    #push!(selectedcols, (287, "K"))
    #push!(selectedcols, (152, "P"))
    push!(selectedcols, (438, "I"))
    #push!(selectedcols, (198, "T"))
    #push!(selectedcols, (214, "G"))
    #push!(selectedcols, (502, "A"))
    #push!(selectedcols, (510, "Q"))
    push!(selectedcols, (384, "N"))

    singlecolumns = false
    if singlecolumns
        fout = open("substitutionmapping.csv","w")
        println(fout,"Mutation,freq(AA),freq(AA|LP),freq(AA|HP),freq(LP|AA),freq(HP|AA),entropy,discrimination")
        for x=1:numcols
            col1 = x
            #col1 = selectedcols[x][1]
            #col1aa = selectedcols[x][2]
            aatime = zeros(Float64,2,20)
            for z=1:numsamples
                allpaths = Array{Array{Int,1},1}[]
                alltimes = Array{Array{Float64,1},1}[]
                selcols = Int[col1]
                for col in selcols
                    nodepaths = Array{Int,1}[]
                    nodetimes = Array{Float64,1}[]
                    push!(allpaths, nodepaths)
                    push!(alltimes, nodetimes)
                    index = CommonUtils.sample(rng, sitelikelihoodsconditionals[col,:])
                    logm = zeros(Float64, length(nodelist))
                    aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)
                    Q = minx[1]*LGmatrix*rates[index]
                    sample = backwardssampling(rng, nodelist, aalikelihoods, Q, LGfreqs, 1)
                    push!(nodepaths,Int[])
                    push!(nodetimes,Float64[])
                    for node in nodelist
                        if !isroot(node)
                            parent = get(node.parent)
                            path,time = CTMCs.modifiedrejectionsampling(rng, Q*node.branchlength, sample[parent.nodeindex], sample[node.nodeindex], 0)
                            push!(nodepaths,path)
                            push!(nodetimes,time)
                        end
                    end
                end

                nodepaths = Array{Int,1}[]
                nodetimes = Array{Float64,1}[]
                push!(allpaths, nodepaths)
                push!(alltimes, nodetimes)
                T,freqs = getTmatrix(initial.tau, Float64[initial.p, 1.0-initial.p])
                traitsample = backwardssampling(rng, nodelist, traitlikelihoods, T, freqs, 1)
                push!(nodepaths,Int[])
                push!(nodetimes,Float64[])
                for node in nodelist
                    if !isroot(node)
                        parent = get(node.parent)
                        path,time = CTMCs.modifiedrejectionsampling(rng, T*node.branchlength, traitsample[parent.nodeindex], traitsample[node.nodeindex], 0)
                        push!(nodepaths,path)
                        push!(nodetimes,time)
                    end
                end

                branchpaths = BranchPath[]
                for node in nodelist
                    colpaths = Array{Int,1}[]
                    coltimes = Array{Float64,1}[]
                    for col=1:2
                        push!(colpaths, allpaths[col][node.nodeindex])
                        push!(coltimes, alltimes[col][node.nodeindex])
                    end
                    branchpath = BranchPath(colpaths, coltimes)
                    push!(branchpaths, branchpath)
                end

                for node in nodelist
                    if !isroot(node)
                        branchiterator = BranchPaths.BranchPathIterator(branchpaths[node.nodeindex],Int[1,2])
                        for (prevstates, prevtime, currstates, currtime, mincol) in branchiterator
                            aa1 = prevstates[1]
                            trait = prevstates[2]
                            time = node.branchlength*(currtime-prevtime)
                            aatime[trait,aa1] += time
                        end
                    end
                end
            end


            aafreqs = sum(aatime,dims=1)/sum(aatime)
            aatimeLP = aatime[1,:]/sum(aatime[1,:])
            aatimeHP = aatime[2,:]/sum(aatime[2,:])
            empiricaltraitfreqs = vec(sum(aatime,dims=2)/sum(aatime))
            #println(col1,"\t", aminoacids[aa1index],"\t", @sprintf("%0.5f", empiricaltraitfreqs[1]), "\t", @sprintf("%0.5f", empiricaltraitfreqs[2]))
            #println("AA\tfreqAA\tAA|LP\tAA|HP\tLP|AA\tHP|AA")

            discrimination = 0.0
            for aa=1:20
                aatimeaasum = sum(aatime[:,aa])
                HPcondAA = 0.0
                if aafreqs[aa] > 0.0
                    HPcondAA = aatime[2,aa]/aatimeaasum
                end
                discrimination += aafreqs[aa]*max(HPcondAA,1.0-HPcondAA)
            end
            informationentropy = 0.0
            for x=1:size(aatime,1)
                for y=1:size(aatime,2)
                    prob = aatime[x,y]/sum(aatime)
                    if prob > 0.0
                        informationentropy -= prob*log(prob)
                    end
                end
            end
            aaindex = 1
            for (all,LP,HP) in zip(aafreqs,aatimeLP,aatimeHP)
                aatimeaasum = sum(aatime[:,aaindex])
                aatimeaa = vec(aatime[:,aaindex]/aatimeaasum)
                if aatimeaasum == 0.0
                    aatimeaa = zeros(2)
                end

                if all > 0.0
                    println(fout,"$(col1)$(aminoacids[aaindex]),", @sprintf("%0.5f", all), ",", @sprintf("%0.5f", LP), ",", @sprintf("%0.5f", HP), ",", @sprintf("%0.5f", aatimeaa[1]), ",", @sprintf("%0.5f", aatimeaa[2]),",",informationentropy,",",discrimination)
                    #println("$(col1)$(aminoacids[aaindex])\t", @sprintf("%0.5f", HP), "\t", @sprintf("%0.5f", aatimeaa[2]),"\t",informationentropy)
                    flush(fout)
                end
                aaindex += 1
            end
        end
        close(fout)
    end

    for x=1:length(selectedcols)
        col1 = selectedcols[x][1]
        col1aa = selectedcols[x][2]
        for y=x+1:length(selectedcols)
            col2 = selectedcols[y][1]
            col2aa = selectedcols[y][2]
            aa1index = findfirst(col1aa, aminoacids)[1]
            aa2index = findfirst(col2aa, aminoacids)[1]

            coevolutionmatrixtotal = zeros(2,20,20)
            observedvals = Float64[]
            expectedvals = Float64[]
            for z=1:numsamples
                coevolutionmatrixsample = zeros(2,20,20)
                allpaths = Array{Array{Int,1},1}[]
                alltimes = Array{Array{Float64,1},1}[]
                selcols = Int[col1,col2]
                for col in selcols
                    nodepaths = Array{Int,1}[]
                    nodetimes = Array{Float64,1}[]
                    push!(allpaths, nodepaths)
                    push!(alltimes, nodetimes)
                    index = CommonUtils.sample(rng, sitelikelihoodsconditionals[col,:])
                    logm = zeros(Float64, length(nodelist))
                    aalikelihoods = getaalikelihoods(nodelist, sequences, getaacolumn(sequences, col), seqindextonodeindex)
                    Q = minx[1]*LGmatrix*rates[index]
                    sample = backwardssampling(rng, nodelist, aalikelihoods, Q, LGfreqs, 1)
                    push!(nodepaths,Int[])
                    push!(nodetimes,Float64[])
                    for node in nodelist
                        if !isroot(node)
                            parent = get(node.parent)
                            path,time = CTMCs.modifiedrejectionsampling(rng, Q*node.branchlength, sample[parent.nodeindex], sample[node.nodeindex], 0)
                            push!(nodepaths,path)
                            push!(nodetimes,time)
                        end
                    end
                end

                nodepaths = Array{Int,1}[]
                nodetimes = Array{Float64,1}[]
                push!(allpaths, nodepaths)
                push!(alltimes, nodetimes)
                T,freqs = getTmatrix(initial.tau, Float64[initial.p, 1.0-initial.p])
                traitsample = backwardssampling(rng, nodelist, traitlikelihoods, T, freqs, 1)
                push!(nodepaths,Int[])
                push!(nodetimes,Float64[])
                for node in nodelist
                    if !isroot(node)
                        parent = get(node.parent)
                        path,time = CTMCs.modifiedrejectionsampling(rng, T*node.branchlength, traitsample[parent.nodeindex], traitsample[node.nodeindex], 0)
                        push!(nodepaths,path)
                        push!(nodetimes,time)
                    end
                end

                branchpaths = BranchPath[]
                for node in nodelist
                    colpaths = Array{Int,1}[]
                    coltimes = Array{Float64,1}[]
                    for col=1:3
                        push!(colpaths, allpaths[col][node.nodeindex])
                        push!(coltimes, alltimes[col][node.nodeindex])
                    end
                    branchpath = BranchPath(colpaths, coltimes)
                    push!(branchpaths, branchpath)
                end

                for node in nodelist
                    if !isroot(node)
                        branchiterator = BranchPaths.BranchPathIterator(branchpaths[node.nodeindex],Int[1,2,3])
                        for (prevstates, prevtime, currstates, currtime, mincol) in branchiterator
                            aa1 = prevstates[1]
                            aa2 = prevstates[2]
                            trait = prevstates[3]
                            time = node.branchlength*(currtime-prevtime)
                            coevolutionmatrixtotal[trait,aa1,aa2] += time
                            coevolutionmatrixsample[trait,aa1,aa2] += time
                        end
                    end
                end
                coevolutionmatrixsample = coevolutionmatrixsample[2,:,:]
                coevolutionmatrixsample /= sum(coevolutionmatrixsample)
                marginalaa1sample = vec(sum(coevolutionmatrixsample,dims=2))
                marginalaa2sample = vec(sum(coevolutionmatrixsample,dims=1))
                observedval = coevolutionmatrixsample[aa1index,aa2index]
                expectedval = marginalaa1sample[aa1index]*marginalaa2sample[aa2index]
                #println(observedval,"\t",expectedval)
                push!(observedvals, observedval)
                push!(expectedvals, expectedval)
            end
            difference = observedvals./expectedvals
            #println(observedvals./expectedvals)
            println("median ", median(difference))
            println("lower ", quantile(difference,0.025))
            println("upper ", quantile(difference,0.975))


            coevolutionobservedall = sum(coevolutionmatrixtotal, dims=1)[1,:,:]
            coevolutionobservedall /= sum(coevolutionobservedall)
            marginalaa1all = vec(sum(coevolutionobservedall,dims=2))
            marginalaa2all = vec(sum(coevolutionobservedall,dims=1))
            coevolutionexpected = zeros(Float64, 20, 20)
            for x=1:20
                for y=1:20
                    coevolutionexpected[x,y] = marginalaa1all[x]*marginalaa2all[y]
                end
            end
            println("All $(col1) $(col1aa) and $(col2) $(col2aa)\t",coevolutionobservedall[aa1index,aa2index],"\t", coevolutionexpected[aa1index,aa2index],"\t", marginalaa1all[aa1index],"\t", marginalaa2all[aa2index])

            coevolutionobservedLP = coevolutionmatrixtotal[1,:,:]
            coevolutionobservedLP /= sum(coevolutionobservedLP)
            marginalaa1LP = vec(sum(coevolutionobservedLP,dims=2))
            marginalaa2LP = vec(sum(coevolutionobservedLP,dims=1))
            coevolutionexpected = zeros(Float64, 20, 20)
            for x=1:20
                for y=1:20
                    coevolutionexpected[x,y] = marginalaa1LP[x]*marginalaa2LP[y]
                end
            end
            println("LP $(col1) $(col1aa) and $(col2) $(col2aa)\t",coevolutionobservedLP[aa1index,aa2index],"\t", coevolutionexpected[aa1index,aa2index],"\t", marginalaa1LP[aa1index],"\t", marginalaa2LP[aa2index])

            coevolutionobservedHP = coevolutionmatrixtotal[2,:,:]
            coevolutionobservedHP /= sum(coevolutionobservedHP)
            marginalaa1HP = vec(sum(coevolutionobservedHP,dims=2))
            marginalaa2HP = vec(sum(coevolutionobservedHP,dims=1))
            coevolutionexpected = zeros(Float64, 20, 20)
            for x=1:20
                for y=1:20
                    coevolutionexpected[x,y] = marginalaa1HP[x]*marginalaa2HP[y]
                end
            end
            println("HP $(col1) $(col1aa) and $(col2) $(col2aa)\t",coevolutionobservedHP[aa1index,aa2index],"\t", coevolutionexpected[aa1index,aa2index],"\t", marginalaa1HP[aa1index],"\t", marginalaa2HP[aa2index])
            println("------------------------------------------------")
        end
    end
end

function loadtree(treefile)
    newickin = open(treefile,"r")
    newickstring = strip(readlines(newickin)[1])
    close(newickin)
    root = gettreefromnewick(newickstring)
    doroottree = true
    if length(root.children) == 2
        doroottree = false
    end
    if doroottree
        root = roottree(gettreefromnewick(newickstring), 1)
    end
    nodelist = getnodelist(root)
end

function simulation()
    parsed_args = parse_commandline()
    collapsetraits = Int[]
    if parsed_args["collapsetraits"] != nothing
        collapsetraits = Int[parse(Int,s) for s in split(parsed_args["collapsetraits"],",")]
    end
    return simulation(parsed_args["alignment"], parsed_args["tree"], split(parsed_args["annotations"],","), collapsetraits, parsed_args["output"])
end

function processsimulations()
    simrates = Float64[0.4,1.0,2.5]

    #simrates = Float64[0.4]
    taus = [2.0, 4.0, 7.5]
    #taus = [2.0]
    lambdas = [2.0,4.0, 8.0, 16.0]
    #lambdas = [2.0,4.0]
    sims = 1:7

    for mu in simrates
        for tau in taus
            for lambda in lambdas
                P_ml = 0.0
                N_ml = 0.0
                TP_ml = 0.0
                FN_ml = 0.0
                FP_ml = 0.0
                TN_ml = 0.0
                total_ml = 0.0

                P_bh1 = 0.0
                N_bh1 = 0.0
                TP_bh1 = 0.0
                FN_bh1 = 0.0
                FP_bh1 = 0.0
                TN_bh1 = 0.0

                for sim in sims
                    infile = "simfinal.fas.lambda$(lambda).rate$(mu).tau$(tau).sim$(sim).fas.ml.csv"
                    fin = open(infile, "r")

                    lineno = 0
                    for line in readlines(fin)
                        spl = split(line,",")
                        if lineno > 0
                            test = div(lineno-1,21) + 1
                            aa = (lineno-1) % 21

                            chi2_ml = parse(Float64, spl[11])
                            pval_ml = parse(Float64, spl[12])
                            #=
                            if chi2_ml > 0.0 && test == aa
                                if pval_ml < 0.05
                                    TP_ml += 1.0
                                else
                                    FN_ml += 1.0
                                end
                                P_ml += 1.0
                            elseif chi2_ml > 0.0
                                if pval_ml < 0.05
                                    FP_ml += 1.0
                                else
                                    TN_ml += 1.0
                                end
                                N_ml += 1.0
                            end
                            =#
                            if chi2_ml > 0.0 && test == aa
                                if pval_ml < 0.05
                                    TP_ml += 1.0
                                else
                                    FN_ml += 1.0
                                end
                        elseif aa > 0 && test  > 20
                                    if pval_ml< 0.05 && chi2_ml > 0.0
                                        FP_ml += 1.0
                                    else
                                        TN_ml += 1.0
                                    end
                            end


                            chi2_bh1 = parse(Float64, spl[17])
                            pval_bh1 = parse(Float64, spl[18])
                            if chi2_bh1 > 0.0 && test == aa
                                if pval_bh1 < 0.05
                                    TP_bh1 += 1.0
                                else
                                    FN_bh1 += 1.0
                                end
                           elseif aa > 0 && test  > 20
                                    if pval_bh1 < 0.05 && chi2_bh1 > 0.0
                                        FP_bh1 += 1.0
                                    else
                                        TN_bh1 += 1.0
                                    end
                            end
                            #=
                            elseif chi2_bh1 > 0.0
                                if pval_bh1 < 0.05
                                    FP_bh1 += 1.0
                                else
                                    TN_bh1 += 1.0
                                end
                                N_bh1 += 1.0
                            end=#
                        end
                        lineno += 1
                    end
                end

                TPR_ml = TP_ml/(TP_ml+FN_ml)
                recall_ml = TP_ml/(TP_ml+FN_ml+1e-10)
                precision_ml = TP_ml/(TP_ml+FP_ml+1e-10)
                TNR_ml = TN_ml/(TN_ml+FP_ml)

                TPR_bh1 = TP_bh1/(TP_bh1+FN_bh1)
                recall_bh1 = TP_bh1/(TP_bh1+FN_bh1+1e-10)
                precision_bh1 = TP_bh1/(TP_bh1+FP_bh1+1e-10)
                TNR_bh1 = TN_bh1/(TN_bh1+FP_bh1)
                #println("TN ", TN_ml,"\t",TN_bh1,"\t",FP_ml,"\t",FP_bh1)
                #println(mu,"\t",tau,"\t",lambda,"\t",@sprintf("%0.4f", TPR_ml),"\t",@sprintf("%0.4f", TNR_ml)),"\t",@sprintf("%0.4f", TPR_bh1),"\t",@sprintf("%0.4f", TNR_bh1))

                println(lambda,"\t",tau,"\t",@sprintf("%0.2f", recall_ml),"\t",@sprintf("%0.2f", precision_ml),"\t",@sprintf("%0.2f", fscore(recall_ml,precision_ml)),"\t",@sprintf("%0.2f", recall_bh1),"\t",@sprintf("%0.2f", precision_bh1),"\t",@sprintf("%0.2f", fscore(recall_bh1,precision_bh1)))
            end
        end
    end

end

function simulation(fastafile::AbstractString, treefile, annotationnames, collapsetraits::Array{Int,1}, outfilein::AbstractString)
    if treefile == nothing
        treestring, treefile = Binaries.fasttreeaa(fastafile)
    end

    rng = MersenneTwister(2184104820249809)
    zup = 1

    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
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

    nodelist = loadtree(treefile)
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
    println(seqnametonodeindex)
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

    #collapsetraits = Int[1,2]
    nodetraits = zeros(Int, length(nodelist))
    for i=1:length(traits)
        if !(traits[i] in collapsetraits)
            nodetraits[seqindextonodeindex[i]] = traits[i]
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
            traitvalue =  argmax(traitcounts[nodeindex,:])
            markednodes[nodeindex] = traitvalue
            if traitvalue in collapsetraits
                nodetraits[nodeindex] = traitvalue
            end
        end
    end
    #println("markednodes ", markednodes)
    #println("nodetraits  ", nodetraits)

    initial = AATraitParameters(Float64[31.47925915,1.0,7.456053451, 0.394240064])
    traitlikelihoods = gettraitlikelihoods(nodelist, traits, seqindextonodeindex)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_DIRECT)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :GN_CRS2_LM)
    println(getarray(initial))
    initial = optimizetraitmodel(nodelist, traitlikelihoods, initial, :LN_COBYLA)
    println(getarray(initial))

    numsamples = 500
    numcats = 10
    aadata,subcolumnrefs = getaadata(nodelist,seqindextonodeindex,sequences)
    initialmu = 0.834965
    initialshape =  0.0936208
    minf, minx = optimizealignmentlikelihood(nodelist, numcols, aadata, subcolumnrefs,initialmu,initialshape,numcats, 100)
    println("Finished",minx)
    initial.mu =  minx[1]
    println("MU ",minx[1])
    maxll, sitelikelihoodsconditionals = aa_alignmentloglikelihood(nodelist, numcols, aadata, subcolumnrefs, minx[1], minx[2], numcats)
    rates = discretizegamma(minx[2], 1.0/minx[2], numcats)
    println("RATES",rates)
    startmu = initial.mu

    simrates = Float64[0.4,1.0,2.5]
    #numcats = length(rates)

    numseqs = 0
    for node in nodelist
        if isleafnode(node)
            numseqs +=1
        end
    end

    traitdata = zeros(Int, numseqs)
    T,freqs = getTmatrix(initial.tau, Float64[initial.p, 1.0-initial.p])
    traitlikelihoods = ones(length(nodelist),2)
    sample = backwardssampling(rng, nodelist, traitlikelihoods, T, freqs, 1)
    for node in nodelist
        if isleafnode(node)
            #traitdata[node.seqindex] = sample[node.nodeindex]
        end
    end
    #traitdata = copy(collapsetraits)
    traitdata = copy(traits)

    aadata = zeros(Int, numseqs, numsamples)
    for z=1:numsamples
        aalikelihoods = ones(length(nodelist),20)
        Q = minx[1]*LGmatrix*rand(rng,rates)
        sample = backwardssampling(rng, nodelist, aalikelihoods, Q, LGfreqs, 1)

        for node in nodelist
            if isleafnode(node)
                aadata[node.seqindex,z] = sample[1,node.nodeindex]
            end
        end
    end

    #=
    aadata = zeros(Int, numseqs, numsamples)
    twocounts = 0.0
    twototals = 0.0
    for z=1:numsamples
        counts = Dict{Int,Int}()
        aalikelihoods = ones(length(nodelist),20)
        z = rand(rng, 1:length(rates))
        Q = rates[z]*minx[1]*LGmatrix*rand(rng,rates)
        sample = backwardssampling(rng, nodelist, aalikelihoods, Q, LGfreqs, 1)

        for node in nodelist
            if isleafnode(node)
                aadata[node.seqindex,z] = sample[1,node.nodeindex]
                counts[sample[1,node.nodeindex]] = get(counts, sample[1,node.nodeindex], 0) + 1
            end
        end

        if length(counts) > 1
            twocounts += 1.0
        end
        twototals += 1.0

        println(length(counts),"\t",counts,"\t", twocounts/twototals)
    end

    exit()=#

    independentsims = false
    taus = [2.0, 4.0, 7.5]
    lambdas = [2.0,4.0, 8.0, 16.0]
    for simulationno=1:1000
        for tau in taus
            for lambda in lambdas
                for catno=1:length(simrates)
                    aatraitdata = zeros(Int, numseqs, numsamples)
                    aarate = simrates[catno]
                    outfile = string(outfilein,".lambda",lambda,".rate",aarate,".tau",tau)
                    for z=1:numsamples
                        targetaa = (z-1) % 20 + 1

                        initial.tau = 1.0
                        if z <= 20
                            initial.mu = tau
                            initial.lambda = lambda
                            aarate = simrates[catno]
                        else
                            initial.mu = startmu
                            initial.lambda = 1.0
                            aarate = rates[(z-1)%numcats+1]
                        end
                        #aarate = 0.1
                        println(z,"\t",initial.lambda,"\t",aarate,"\t",targetaa)
                        Q,freqs,marginaltraitfreqs = getAATraitmatrix(initial, LGmatrix*aarate, LGfreqs, targetaa)

                        #aatraitlikelihoods = getaatraitlikelihoods(nodelist, sequences, zeros(Int,length(sequences)), traitdata, seqindextonodeindex)
                        aatraitlikelihoods = zeros(length(nodelist),40)
                        for seqindex=1:numseqs
                            for aaindex=1:20
                                aatraitlikelihoods[seqindextonodeindex[seqindex], (traitdata[seqindex]-1)*20 + aaindex] = 1.0
                            end
                        end
                        sample = backwardssampling(rng, nodelist, aatraitlikelihoods, Q, freqs, 1)
                        for node in nodelist
                            if isleafnode(node)
                                aatraitdata[node.seqindex,z] = sample[node.nodeindex]
                                trait = div(aatraitdata[node.seqindex,z]-1,20) + 1
                                sampledaa = ((aatraitdata[node.seqindex,z]-1) % 20) + 1
                                aadata[node.seqindex,z] = sampledaa
                            end
                        end
                    end

                    outsequences = AbstractString[]
                    outnames = AbstractString[]
                    for seqindex=1:numseqs
                        seq = ""
                        for z=1:numsamples
                            seq = string(seq, aminoacids[aadata[seqindex,z]])
                        end
                        name = ""
                        if traitdata[seqindex] == 2
                            name = "seq$(seqindex)_HP"
                        elseif traitdata[seqindex] == 1
                            name = "seq$(seqindex)_LP"
                        end
                        nodeindex = seqindextonodeindex[seqindex]
                        nodelist[nodeindex].name = name
                        push!(outnames, name)
                        push!(outsequences, seq)
                    end

                    fastafile = string(outfile,".sim",simulationno,".fas")
                    fout = open(fastafile,"w")
                    for seqindex=1:numseqs
                        write(fout,">$(outnames[seqindex])\n")
                        write(fout,"$(outsequences[seqindex])\n")
                    end
                    close(fout)

                    treefile = string(fastafile,".nwk")
                    fout = open(treefile,"w")
                    println(fout,getnewick(nodelist[1]))
                    close(fout)



                    annotations = "LP,HP"
                    maxcols = 40

                    #=
                    outfilename = string(fastafile,".csv")
                    mlmodel(rng, fastafile,treefile,outfilename,annotations,maxcols)
                    =#

                    mlnewick,cachefile = Binaries.fasttreeaa(fastafile)
                    mlnewickfile = string(treefile,".ml.nwk")
                    fout = open(mlnewickfile,"w")
                    write(fout, mlnewick)
                    close(fout)

                    outfilename = string(fastafile,".ml.csv")
                    mlmodel(rng, fastafile,mlnewickfile,outfilename,annotations,maxcols)
                end
            end
        end
    end
    #=
    T,freqs = getTmatrix(initial.tau, Float64[initial.p, 1.0-initial.p])
    traitsample = backwardssampling(rng, nodelist, traitlikelihoods, T, freqs, 1)
    =#



    #println(outsequences)

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

function marginaltraitprobability(aa, targetaa, piHP::Float64, lambda::Float64, freqs::Array{Float64,1}=LGfreqs)
    #marginalHP = 
    marginalHP = piHP*(lambda + (1.0/lambda))/((piHP*(lambda + (1.0/lambda))) + 2.0*(1.0-piHP))
    #marginalHP = piHP*(lambda + (1.0/lambda))/((piHP*(lambda + (1.0/lambda))) + 2.0*(1.0-piHP))
    #println(aa,"\t",targetaa,"\t", piHP, "\t", lambda,"\t",marginalHP)
    if aa == targetaa
        return piHP*lambda/((1.0-piHP) + piHP*lambda)
    else
        return piHP/((1.0-piHP)*lambda + piHP)
    end
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

function getcolumn(csvfile, key::String)
    fin = open(csvfile, "r")
    header = ""
    data = []
    colindex = 1
    for line in readlines(fin)
        spl = split(strip(string(line)), [',', '\"'],keepempty=false)
        if header == ""
            header = spl
            for (i,colname) in  enumerate(spl)
                #println("!",strip(colname),",",strip(key))
                if colname == key
                    colindex = i
                end
            end
        else
            push!(data, spl[colindex])
        end
    end
    close(fin)
    integer = true
    numeric = true
    for v in data
        try
            parse(Int,v)
        catch
            integer = false
        end
        try
            parse(Float64,v)
        catch
            numeric = false
        end
    end
    if integer
        retdata = Int[]
        for v in data
            push!(retdata, parse(Int,v))
        end
        return retdata
    elseif numeric
        retdata = Float64[]
        for v in data
            push!(retdata, parse(Float64,v))
        end
        return retdata
    end
    return data

end

include("Mapping.jl")
function scoresequence2(mappingsequence::String, segment="H7HA")
    alignment = ""

    scores = []
    if segment == "H7HA"
        alignment = "../data/H7NX/HA/H7_HA_alnP.fas"
        push!(scores, (143, 'A', 0.36327, 0.47780))
        push!(scores, (143, 'T', 0.28351, 7.04113))
        push!(scores, (274, 'I', 0.16081, 6.19619))
        push!(scores, (438, 'I', 0.17541, 4.62822))
    end
    mapping, revmapping = Mapping.createmapping(alignment, mappingsequence)
    v = Float64[0.5, 0.5]
    for (pos,targetaa,p,lambda) in scores
        seqaa = get(mappingsequence, mapping[pos], "-")
        marginal = marginaltraitprobability(seqaa, targetaa, p, lambda)
        v[1] *= marginal
        v[2] *= (1.0-marginal)
        v /= sum(v)
        println(pos,"\t",seqaa,"\t",targetaa,"\t",v,"\t",marginal)
    end
end

function geommean(probs::Array{Float64,1})
    wi = 1.0
    u = sum(log.(probs.^wi))
    v = logsumexp(u, sum((1.0 .- probs).^wi))
    return exp(u - v)
end

function satopaaprobcombinations(probs::Array{Float64,1})
    v = 0.0
    wi = 1.0
    alpha = 1.0
    for p in probs
        v += log((p/(1.0-p))^wi)
    end
    v *= alpha
    return exp(v - logsumexp(0.0, v))
end

#=
function scoresequence(mappingsequence::String, segment="H7HA")
    alignment = ""

    targetpositions = [("HA",143),("HA",274),("HA", 438),("HA",384)]
    segments = AbstractString[]
    positions = Int[]
    targets = AbstractString[]
    ps = Float64[]
    lambdas = Float64[]
    if segment == "H7HA"
        alignment = "../data/H7NX/HA/H7_HA_alnP.fas"        
        csvfile = "../data/H7NX/HA/H7_HA_alnP.fas.results.update - Copy (13).csv"
        segments = ParallelEvolution.getcolumn(csvfile, "Segment")
        positions = ParallelEvolution.getcolumn(csvfile, "Site")
        targets = ParallelEvolution.getcolumn(csvfile, "Target AA")
        ps = ParallelEvolution.getcolumn(csvfile, "p")
        lambdas = ParallelEvolution.getcolumn(csvfile, "lambda")        
    end

    marginalprobs = Float64[]
    mapping, revmapping = Mapping.createmapping(alignment, mappingsequence)
    for (segment,pos,target,p,lambda) in zip(segments,positions,targets,ps,lambdas)
        for (targetsegment, seqpos) in targetpositions
            #targetsegment == segment && 
            if seqpos == pos && mapping[pos] > 0 && mappingsequence[mapping[seqpos]] == target[1]
                marginalprob = marginaltraitprobability(target[1], target[1], p, lambda)
                push!(marginalprobs, marginalprob)
                println(segment,"\t",pos,"\t",target,"\t",p,"\t",lambda, "\tmarginal: ", marginalprob)
            end
        end
    end
    println(marginalprobs)
    println(satopaaprobcombinations(marginalprobs))
    #println(geommean(marginalprobs))
    
    #=
    mapping, revmapping = Mapping.createmapping(alignment, mappingsequence)
    v = Float64[0.5, 0.5]
    for (pos,targetaa,p,lambda) in scores
        seqaa = get(mappingsequence, mapping[pos], "-")
        marginal = marginaltraitprobability(seqaa, targetaa, p, lambda)
        v[1] *= marginal
        v[2] *= (1.0-marginal)
        v /= sum(v)
        println(pos,"\t",seqaa,"\t",targetaa,"\t",v,"\t",marginal)
    end=#
end

function scorealignment(fastafile, segment::AbstractString="H7HA")
    sequences = AbstractString[]
    names = AbstractString[]
    FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names, replace(strip(desc), "\\r" => ""))
            push!(sequences, replace(seq, "\r" => ""))
        end
    end
    for (name,seq) in zip(names,sequences)
        println(name)
        scoresequence(seq, segment)
    end

end
=#

function scoresequence(alignmentfile::AbstractString, mappingsequence::AbstractString, csvfile::AbstractString, targetpositions::Array{Int,1}, mapping=nothing, revmapping=nothing)
    segments = AbstractString[]
    positions = Int[]
    targets = AbstractString[]
    ps = Float64[]
    lambdas = Float64[]

    segments = ParallelEvolution.getcolumn(csvfile, "Segment")
    positions = ParallelEvolution.getcolumn(csvfile, "Site")
    targets = ParallelEvolution.getcolumn(csvfile, "Target AA")
    ps = ParallelEvolution.getcolumn(csvfile, "p")
    lambdas = ParallelEvolution.getcolumn(csvfile, "lambda")

    marginalprobs = Float64[]
    if mapping == nothing
        println("CREATING MAPPING")
        mapping, revmapping = Mapping.createmapping(alignmentfile, mappingsequence)
    end
    
    for (segment,pos,target,p,lambda) in zip(segments,positions,targets,ps,lambdas)
        for seqpos in targetpositions
            #targetsegment == segment && 
            if seqpos == pos && mapping[pos] > 0 && mappingsequence[mapping[seqpos]] == target[1]
                marginalprob = marginaltraitprobability(target[1], target[1], p, lambda)
                push!(marginalprobs, marginalprob)
                #println(segment,"\t",pos,"\t",target,"\t",p,"\t",lambda, "\tmarginal: ", marginalprob)
            end
        end
    end

    #=
    targetpositions = [(143,'T'),(274,'I'),(438,'I'),(384,'N')]
    for (segment,pos,target,p,lambda) in zip(segments,positions,targets,ps,lambdas)
        for (seqpos,targetaa) in targetpositions
            #targetsegment == segment && 
            if seqpos == pos && mapping[pos] > 0 && targetaa == target[1] && mappingsequence[mapping[seqpos]] == target[1]
                marginalprob = marginaltraitprobability(target[1], target[1], p, lambda)
                push!(marginalprobs, marginalprob)
                #println(segment,"\t",pos,"\t",target,"\t",p,"\t",lambda, "\tmarginal: ", marginalprob)
            end
        end
    end=#
    return marginalprobs
end

function collatesequences(alignments::Array{AbstractString,1}, segmentnames::Array{AbstractString,1})
    seqdict = Dict{String, Array{Tuple{AbstractString,AbstractString},1}}()
    for (segment,alignment) in zip(segmentnames,alignments)
        if isfile(alignment)
            sequences = AbstractString[]
            names = AbstractString[]
            FastaIO.FastaReader(alignment) do fr
                for (desc, seq) in fr
                    push!(names, replace(strip(desc), "\\r" => ""))
                    push!(sequences, replace(seq, "\r" => ""))
                end
            end
            for (name,seq) in zip(names,sequences)
                canonicalname = replace(name, "|" => ".")
                ls = get(seqdict,canonicalname, Tuple{AbstractString,AbstractString}[])
                push!(ls, (segment,seq))
                seqdict[canonicalname] = ls
            end
        end
    end
    return seqdict
end

function scorelisttomatrix(ls::Array{Array{Float64,1},1})
    mat = zeros(Float64,length(ls),length(ls[1]))
    for (i,arr) in enumerate(ls)
        for (j,v) in enumerate(arr)
            mat[i,j] = v
        end
    end
    return mat
end


using Plots
#using ROC
function main2()
    alignments = AbstractString["../data/H7NX/HA/H7_HA_alnP.fas", "../data/H7NX/PB1/H7_PB1_alnP.fasta", "../data/H7NX/PB2/H7_PB2_alnP.fasta", "../data/H7NX/M/H7_M1_alnP.fasta", "../data/H7NX/NP/H7_NP_alnP.fasta", "../data/H7NX/NS/H7_NS1_alnP.fasta"]
    segments = AbstractString["H7HA","H7PB1", "H7PB2", "H7M", "H7NP", "H7NS1"]
    alignmentdict = Dict{String, String}()
    for (segmentname,alignmentfile) in zip(segments,alignments)
        alignmentdict[segmentname] = alignmentfile
    end
    #lignmentdict = Dict{String, String}()
    #lignmentdict["H7HA"] = "../data/H7NX/HA/H7_HA_alnP.fas"
    #lignmentdict["H7HA"] = "../data/H7NX/HA/H7_HA_alnP.fas"

    csvdict = Dict{String, String}()
    csvdict["H7HA"] = "../data/H7NX/HA/H7_HA_alnP.fas.results.update - Copy (14).csv"
    csvdict["H7PB2"] = "../data/H7NX/PB2/H7_PB2_alnP.fasta.results.update - Copy.csv"
    csvdict["H7M"] = "../data/H7NX/M/H7_M1_alnP.fasta.results.update.csv"
    csvdict["H7NP"] = "../data/H7NX/NP/H7_NP_alnP.fasta.results.update.csv"
    posdict = Dict{String, Array{Int,1}}()
    posdict["H7HA"] = Int[143, 274, 438, 384]
    posdict["H7PB2"] = Int[355, 480]
    #posdict["H7PB1"] = Int[154, 152]
    #posdict["H7NS1"] = Int[209]

    mappings = Dict{AbstractString, Tuple{Dict{Int,Int},Dict{Int,Int}}}()
    for segmentname in segments
        if haskey(alignmentdict, segmentname)
            mapping, revmapping = Mapping.createmapping2(alignmentdict[segmentname], alignmentdict[segmentname])
            mappings[segmentname] = (mapping, revmapping)
        end
    end


    seqdict = ParallelEvolution.collatesequences(alignments,segments)
    satoscores = Float64[]
    geomscores = Float64[]
    scorematrix = Array{Float64,1}[]
    labels = Int[]
    seqindex = 0
    for seqname in keys(seqdict)
        allmarginalprobs = Float64[]
        for (segmentname, seq) in seqdict[seqname]
            mapping,revmapping = mappings[segmentname]
            if haskey(posdict, segmentname)
                marginalprobs = ParallelEvolution.scoresequence(alignmentdict[segmentname], seq, csvdict[segmentname], posdict[segmentname], mapping, revmapping)
                #println(segmentname,"\t",marginalprobs)
                append!(allmarginalprobs, marginalprobs)
            end
        end
        while length(allmarginalprobs) < 6
            push!(allmarginalprobs, 0.2)
        end
        if length(allmarginalprobs) > 0
            classlabel = 0
            if occursin(".HP.", seqname)
                classlabel = 1
            end

            push!(scorematrix,allmarginalprobs)
            satoscore = satopaaprobcombinations(allmarginalprobs)
            push!(satoscores, satoscore)
            geomscore = mean(allmarginalprobs)
            push!(geomscores, geomscore)
            
            
            push!(labels, classlabel)
            println(seqname,"\t",classlabel,"\t", satoscore,"\t",geomscore,"\t", length(seqdict[seqname]),"\t", allmarginalprobs)
            satorocdata = roc(satoscores,labels,1)
            geomrocdata = roc(geomscores,labels,1)
            println(AUC(satorocdata),"\t",AUC(geomrocdata))
            seqindex += 1
            if seqindex % 10 == 0
                plot(satorocdata);
                savefig("sato_roc.png")
                plot(geomrocdata);
                savefig("geom_roc.png")

                #println(scorelisttomatrix(scorematrix))
            end
        end
    end
end

using NPZ
function creatematrix(subtype::AbstractString="H5", queryalignments=AbstractString[])
    alignments = AbstractString[]
    segments = AbstractString[]
    alignmentdict = Dict{String, String}()
    posdict = Dict{String, Array{Tuple{Int,Char},1}}()

    if subtype == "H7"
        alignments = AbstractString["../data/H7NX/HA/H7_HA_alnP.fas", "../data/H7NX/PB1/H7_PB1_alnP.fasta", "../data/H7NX/PB2/H7_PB2_alnP.fasta", "../data/H7NX/M/H7_M1_alnP.fasta", "../data/H7NX/NP/H7_NP_alnP.fasta", "../data/H7NX/NS/H7_NS1_alnP.fasta"]
        segments = AbstractString["H7HA","H7PB1", "H7PB2", "H7M", "H7NP", "H7NS1"]    
        posdict["H7HA"] = [(143,'T'),(274,'I'),(438,'I'),(384,'N')]
        posdict["H7PB2"] = [(355,'K'),(480,'I')]
        posdict["H7NS1"] = [(209,'N')]
        posdict["H7PB1"] = [(154,'D'),(152,'L')]
    else
        alignments = AbstractString[]
        segments = AbstractString["H5HA","H5PB1", "H5PB2", "H5M", "H5NS", "H5PA"] 
        
        
        push!(alignments, "../data/H5_MODEL/H5/H5_aln_shortP.fasta")
        posdict["H5HA"] = [(379,'R'), (242,'I'), (145,'L'), (154,'I'), (154,'L'), (127,'L'), (167,'T'), (204,'I'), (157,'P')]        
        push!(alignments, "../data/H5_MODEL/PB1/PB1_shortP.fasta")
        posdict["H5PB1"] = [(113,'I')]
        push!(alignments, "../data/H5_MODEL/PB2/PB2_shortP.fasta")
        posdict["H5PB2"] = [(674,'T'),(451,'T'),(627,'K'),(508,'Q')]
        push!(alignments, "../data/H5_MODEL/M/M_shortP.fasta")
        posdict["H5M"] = [(101,'K')]
        
        #push!(alignments, "../data/H5NX/NS_shortP.fasta")        
        #push!(alignments, "../data/H5NX/PA_shortP.fasta")

        #=
        push!(alignments, "../data/H5_MODEL/H5/H5_aln_shortP.fasta")
        posdict["H5HA"] = [(379,'R'), (242,'I'), (145,'L')]        
        push!(alignments, "../data/H5_MODEL/PB1/PB1_shortP.fasta")
        posdict["H5PB1"] = [(113,'I')]
        push!(alignments, "../data/H5_MODEL/PB2/PB2_shortP.fasta")
        posdict["H5PB2"] = [(674,'T')]
        push!(alignments, "../data/H5_MODEL/M/M_shortP.fasta")
        posdict["H5M"] = [(101,'K')]
        =#

    end
    alltargets = AbstractString[]
    for segmentname in segments
        if haskey(posdict,segmentname)
            for (pos,targetaa) in posdict[segmentname]
                push!(alltargets, string(segmentname[3:end],"_",pos,targetaa))
            end
        end
    end
    println(join(alltargets,"\t"))
    #exit()

    hasquery = true
    if length(queryalignments) == 0
        queryalignments = alignments
        hasquery = false
    end

    mappings = Dict{AbstractString, Tuple{Dict{Int,Int},Dict{Int,Int}}}()
    for (segmentname, alignment, queryalignment) in zip(segments,alignments,queryalignments)
        if isfile(queryalignment)
            println(alignment,"\t",queryalignment,"\tfound:yes")
            mapping, revmapping = Mapping.createmapping2(alignment, queryalignment)
            mappings[segmentname] = (mapping, revmapping)
            println(mappings[segmentname])
        else
            println(alignment,"\t",queryalignment,"\tfound:no")
        end
    end

    seqdict = ParallelEvolution.collatesequences(queryalignments,segments)
    rows = Array{Float64,1}[]
    charrows = AbstractString[]
    labels = Float64[]
    names = AbstractString[]
    for (seqindex, seqname) in enumerate(keys(seqdict))
        hasallsegments = true
        for segmentname in keys(posdict)
            hasmatch = false
            for (seqsegmentname,seq) in seqdict[seqname]
                if seqsegmentname == segmentname
                    hasmatch = true
                    break
                end
            end
            if !hasmatch
                hasallsegments = false
            end
        end
        if true || hasallsegments
            row = Float64[]
            chars = ""
            targetchars = ""
            for segmentname2 in segments
                isfound = false 
                for (segmentname,seq) in seqdict[seqname]
                    if haskey(posdict, segmentname2)
                        mapping,revmapping = mappings[segmentname2]
                        if segmentname == segmentname2
                            for (pos, targetaa) in posdict[segmentname2]
                                #println(segmentname,"\t",pos,"\t", targetaa, "\t", seq[mapping[pos]])
                                if targetaa == seq[mapping[pos]]
                                    push!(row, 1.0)
                                    chars = string(chars, seq[mapping[pos]])
                                    targetchars = string(targetchars,targetaa)
                                else
                                    push!(row, 0.0)
                                    chars = string(chars, seq[mapping[pos]])
                                    targetchars = string(targetchars,targetaa)
                                end
                            end
                            isfound = true
                        end
                    end
                end
                if !isfound
                    if haskey(posdict, segmentname2)
                        for (pos, targetaa) in posdict[segmentname2]
                            push!(row, NaN)
                            chars = string(chars, "?")
                            targetchars = string(targetchars,targetaa)
                        end
                    end
                end
            end
            classlabel = 0
            if occursin(".HP.", seqname)
                push!(labels, 1.0)
                classlabel = 1
            else
                push!(labels, 0.0)
            end
            println(row,"\t",classlabel)
            push!(rows,row)
            #println(chars,"\t",targetchars)
            push!(charrows, chars)
            push!(names, seqname)
        end
    end

    mat = zeros(Float64, length(rows), length(rows[1]))
    for i=1:size(mat,1)
        for j=1:size(mat,2)
            mat[i,j] = rows[i][j]
        end
    end
    namefile = ""
    if hasquery
        namefile = string(subtype,"_query_names.txt")
        npzwrite(string(subtype,"_query_X.npy"), mat)
        npzwrite(string(subtype,"_query_y.npy"), labels)
    else
        namefile = string(subtype,"_training_names.txt")
        npzwrite(string(subtype,"_training_X.npy"), mat)
        npzwrite(string(subtype,"_training_y.npy"), labels)
    end
    fout = open(namefile,"w")
    for (name,chars) in zip(names,charrows)
        println(fout,name,"\t",chars)
    end
    close(fout)
end

function simpleimputation(subtype="H7", targetalignments::Array{AbstractString,1}=AbstractString[])
    alignments = AbstractString[]
    segments = AbstractString[]
    alignmentdict = Dict{String, String}()
    posdict = Dict{String, Array{Tuple{Int,Char},1}}()

    if subtype == "H7"
        #referencealignments = AbstractString["../data/H7NX/HA/H7_HA_alnP.fas", "../data/H7NX/PB1/H7_PB1_alnP.fasta", "../data/H7NX/PB2/H7_PB2_alnP.fasta", "../data/H7NX/M/H7_M1_alnP.fasta", "../data/H7NX/NP/H7_NP_alnP.fasta", "../data/H7NX/NS/H7_NS1_alnP.fasta"]
        referencealignments = AbstractString["../data/H7NX/HA/H7_HA_alnP.fas", "../data/H7NX/PB1/H7_PB1_alnP.fasta", "../data/H7NX/PB2/H7_PB2_alnP.fasta", "../data/H7NX/M/H7_M1_alnP.fasta", "../data/H7NX/NP/H7_NP_alnP.fasta", "../data/H7NX/NS/H7_NS1_alnP.fasta", "", "../data/H7NX/PA/H7_PA_alnP.fasta"]
        segments = AbstractString["H7HA","H7PB1", "H7PB2", "H7M", "H7NP", "H7NS1", "H7NA", "H7PA"]   
        posdict["H7HA"] = [(143,'T'),(274,'I'),(438,'I'),(384,'N')]
        push!(posdict["H7HA"], (402, 'K'))
        push!(posdict["H7HA"], (335, 'T'))
        push!(posdict["H7HA"], (287, 'K'))
        push!(posdict["H7HA"], (175, 'G'))
        push!(posdict["H7HA"], (152, 'P'))
        push!(posdict["H7HA"], (175, 'E'))
      
        posdict["H7PB1"] = [(154,'D'),(152,'L')]
        push!(posdict["H7PB1"], (473, 'I'))
        push!(posdict["H7PB1"], (709, 'I'))

        posdict["H7PB2"] = [(355,'K'),(480,'I')]
        push!(posdict["H7PB2"], (356,'I'))
        push!(posdict["H7PB2"], (584,'I'))
        push!(posdict["H7PB2"], (640,'I'))
        push!(posdict["H7PB2"], (655,'A'))
        push!(posdict["H7PB2"], (661,'A'))
        push!(posdict["H7PB2"], (702,'R'))

        posdict["H7NS1"] = [(209,'N')]

       # posdict["H7M"] = []
        #push!(posdict["H7M"], (95, 'K'), (97, 'I'), (131, 'V'))
        #posdict["H7NP"] = []
        #push!(posdict["H7NP"], (101, 'E'))
        posdict["H7NS1"] = []
        push!(posdict["H7NS1"], (56, 'A'))
        push!(posdict["H7NS1"], (180, 'T'))
        #posdict["H7PA"] = []
        #push!(posdict["H7PA"], (32, 'M'))
        #push!(posdict["H7PA"], (63, 'I'))
        #push!(posdict["H7PA"], (115, 'K'))
        #push!(posdict["H7PA"], (237, 'K'))
    else
        referencealignments = AbstractString[]
        segments = AbstractString["H5HA","H5PB1", "H5PB2", "H5M", "H5NS", "H5PA"] 

        push!(referencealignments, "../data/H5_MODEL/H5/H5_aln_shortP.fasta")
        posdict["H5HA"] = [(379,'R'), (242,'I'), (145,'L'), (154,'I'), (154,'L'), (127,'L'), (167,'T'), (204,'I'), (157,'P')]        
        push!(posdict["H5HA"], (214, 'V'))
        push!(posdict["H5HA"], (199, 'N'))
        push!(posdict["H5HA"], (172, 'T'))
        push!(posdict["H5HA"], (145, 'A'))
        push!(posdict["H5HA"], (145, 'V'))
        push!(posdict["H5HA"], (145, 'I'))
        push!(posdict["H5HA"], (145, 'Q'))
        push!(posdict["H5HA"], (145, 'P'))
        push!(posdict["H5HA"], (519, 'I'))
        push!(posdict["H5HA"], (145, 'D'))
        push!(posdict["H5HA"], (145, 'T'))
        push!(posdict["H5HA"], (170, 'N'))
        push!(posdict["H5HA"], (145, 'N'))        
        push!(posdict["H5HA"], (489,'R'))
        push!(posdict["H5HA"], (87,'V'))
        push!(posdict["H5HA"], (149,'A'))
        push!(posdict["H5HA"], (145, 'E'))  
        #=      
        push!(posdict["H5HA"], (87,'Q'))
        #push!(posdict["H5HA"], (145,'M')) 
        push!(posdict["H5HA"], (200,'I'))
        push!(posdict["H5HA"], (200,'T'))
        push!(posdict["H5HA"], (254,'T'))
        push!(posdict["H5HA"], (127,'L'))
        push!(posdict["H5HA"], (87,'D'))
        push!(posdict["H5HA"], (87,'E'))
        push!(posdict["H5HA"], (200,'P'))
        push!(posdict["H5HA"], (200,'N'))
        push!(posdict["H5HA"], (142,'E'))
        push!(posdict["H5HA"], (488,'I'))
        push!(posdict["H5HA"], (205,'T'))
        push!(posdict["H5HA"], (128,'K'))
        push!(posdict["H5HA"], (200,'E'))
        push!(posdict["H5HA"], (7,'F'))        
        push!(posdict["H5HA"], (205,'E'))        
        push!(posdict["H5HA"], (139,'P'))
        push!(posdict["H5HA"], (87,'M'))
        #push!(posdict["H5HA"], (298,'I'))
        #push!(posdict["H5HA"], (298,'I'))=#
        push!(referencealignments, "../data/H5_MODEL/PB1/PB1_shortP.fasta")
        posdict["H5PB1"] = [(113,'I')]
        #push!(posdict["H5PB1"], (111,'I'))
        push!(referencealignments, "../data/H5_MODEL/PB2/PB2_shortP.fasta")
        posdict["H5PB2"] = [(674,'T'),(451,'T'),(627,'K'),(508,'Q')]
        #push!(posdict["H5PB2"], (494,'I'))
        #push!(posdict["H5PB2"], (649,'I'))
        #push!(posdict["H5PB2"], (495,'I'))
        #push!(posdict["H5PB2"], (108,'S'))
        push!(referencealignments, "../data/H5_MODEL/M/M_shortP.fasta")
        posdict["H5M"] = [(101,'K')]
        push!(posdict["H5M"], (205,'I'))
        #push!(posdict["H5M"], (166,'A'))
    end
    alltargets = AbstractString[]
    for segmentname in segments
        if haskey(posdict,segmentname)
            for (pos,targetaa) in posdict[segmentname]
                push!(alltargets, string(segmentname[3:end],"_",pos,targetaa))
            end
        end
    end
    #println(join(alltargets,"\t"))
    #exit()

    hasquery = true
    if length(targetalignments) == 0
        targetalignments = referencealignments
        hasquery = false
    end

    refseqdict = ParallelEvolution.collatesequences(referencealignments,segments)
    targetseqdict = ParallelEvolution.collatesequences(targetalignments,segments)

    mappings = Dict{AbstractString, Tuple{Dict{Int,Int},Dict{Int,Int}}}()
    for (segmentname, referencealignment, targetalignment) in zip(segments,referencealignments,targetalignments)
        if isfile(targetalignment)
            println(referencealignment,"\t",targetalignment,"\tfound:yes")
            mapping, revmapping = Mapping.createmapping2(referencealignment, targetalignment)
            mappings[segmentname] = (mapping, revmapping)
            println(mappings[segmentname])
        else
            println(referencealignment,"\t",targetalignment,"\tfound:no")
        end
    end

    allseqdict = merge(refseqdict, targetseqdict)


    refseqdict2 = Dict{Tuple{AbstractString,AbstractString}, AbstractString}()
    targetseqdict2 = Dict{Tuple{AbstractString,AbstractString}, AbstractString}()
    for refseqname in keys(refseqdict)     
        for (refsegmentname,refseq) in refseqdict[refseqname]
            refseqdict2[(refseqname,refsegmentname)] = refseq
        end
    end
    for targetseqname in keys(targetseqdict)     
        for (targetsegmentname,targetseq) in targetseqdict[targetseqname]
            targetseqdict2[(targetseqname,targetsegmentname)] = targetseq
        end
    end    
    allseqdict2 = merge(refseqdict2, targetseqdict2)

    for segment in segments
        for targetseqname in keys(targetseqdict)
            if !haskey(targetseqdict2, (targetseqname, segment))            
                bestdist = (0.0,0.0)
                bestseqname = ""
                for refseqname in keys(refseqdict)
                    if haskey(allseqdict2, (refseqname, segment))
                        matches = 0.0
                        total = 1e-100
                        for (refsegment, refseq) in allseqdict[refseqname]
                            if haskey(targetseqdict2, (targetseqname, refsegment))
                                mapping, revmapping = mappings[refsegment]
                                targetseq = targetseqdict2[(targetseqname, refsegment)]
                                for pos=1:length(refseq)
                                    if mapping[pos] > 0
                                        #println(refsegment,"\t",refseq[pos],"\t", targetseq[mapping[pos]])
                                        if refseq[pos] == targetseq[mapping[pos]]
                                            matches += 1.0
                                        end
                                        total += 1.0
                                    end
                                end
                            end
                        end
                        similarity = (matches/total, total)
                        #println(similarity)
                        #println(bestdist)
                        if similarity[1] > bestdist[1] || (similarity[1] == bestdist[1] && similarity[2] > bestdist[2])
                            bestdist = (similarity[1],total)
                            bestseqname = refseqname                    
                           # println(matches,"\t",total,"\t",similarity)
                        end
                    end
                end
                if haskey(mappings, segment)
                    mapping, revmapping = mappings[segment]
                    refseqsegment = allseqdict2[(bestseqname, segment)]
                    reconstructtargetsegment = ""
                    for pos=1:length(revmapping)
                        if revmapping[pos] > 0
                            reconstructtargetsegment = string(reconstructtargetsegment, refseqsegment[revmapping[pos]])
                        else
                            reconstructtargetsegment = string(reconstructtargetsegment, "-")
                        end
                    end
                    targetseqdict2[(targetseqname, segment)] = reconstructtargetsegment
                    println(bestdist)
                    #println("-------------------------------------------------")
                end
            end       
        end
    end

    rows = Array{Float64,1}[]
    charrows = AbstractString[]
    labels = Float64[]
    names = AbstractString[]
    for seqname in keys(targetseqdict)
        row = Float64[]
        chars = ""
        targetchars = ""
        for segmentname2 in segments
            isfound = false
            if haskey(mappings, segmentname2) && haskey(posdict,segmentname2)
                mapping,revmapping = mappings[segmentname2]
                seq = targetseqdict2[(seqname,segmentname2)]
                for (pos, targetaa) in posdict[segmentname2]
                    #println(pos)
                    #println("->", mapping[pos])
                    if mapping[pos] > 0 && targetaa == seq[mapping[pos]]
                        push!(row, 1.0)
                        chars = string(chars, seq[mapping[pos]])
                        targetchars = string(targetchars,targetaa)
                    else
                        push!(row, 0.0)
                        char = 'X'
                        if mapping[pos] > 0
                            char = seq[mapping[pos]]
                        end
                        chars = string(chars, char)
                        targetchars = string(targetchars,targetaa)
                    end
                end
                isfound = true
            end
        end
        classlabel = 0
        println(seqname)
        if occursin(".HP.", seqname)
            push!(labels, 1.0)
            classlabel = 1
        else
            push!(labels, 0.0)
        end
        println(row,"\t",classlabel)
        push!(rows,row)
        #println(chars,"\t",targetchars)
        push!(charrows, chars)
        push!(names, seqname)
    end

    mat = zeros(Float64, length(rows), length(rows[1]))
    for i=1:size(mat,1)
        for j=1:size(mat,2)
            mat[i,j] = rows[i][j]
        end
    end
    namefile = ""
    if hasquery
        namefile = string(subtype,"_query_names.txt")
        npzwrite(string(subtype,"_query_X.npy"), mat)
        npzwrite(string(subtype,"_query_y.npy"), labels)
    else
        namefile = string(subtype,"_training_names.txt")
        npzwrite(string(subtype,"_training_X.npy"), mat)
        npzwrite(string(subtype,"_training_y.npy"), labels)
    end
    fout = open(namefile,"w")
    for (name,chars) in zip(names,charrows)
        println(fout,name,"\t",chars)
    end
    close(fout)

end

end

#AbstractString["H7HA","H7PB1", "H7PB2", "H7M", "H7NP", "H7NS1", "H7NA", "H7PA"]   
#ali

#queryalignments = AbstractString["../data/Metric_Datasets/H7/H7/H7_HA_sampled_P.fasta", "../data/Metric_Datasets/H7/PB1/H7_PB1_sampled_P.fasta", "../data/Metric_Datasets/H7/PB2/H7_PB2_sampled_P.fasta", "../data/Metric_Datasets/H7/M/H7_M_aln_P.fasta", "../data/Metric_Datasets/H7/NP/H7_NP_aln_P.fasta", "../data/Metric_Datasets/H7/NS1/H7_NS1_sampled_P.fasta", "../data/Metric_Datasets/H7/NA/H7_NA_aln_P.fasta", "../data/Metric_Datasets/H7/PA/H7_PA_aln_P.fasta"]
#ParallelEvolution.creatematrix("H7",queryalignments)
#ParallelEvolution.creatematrix("H7")
#ParallelEvolution.simpleimputation("H7")
#ParallelEvolution.simpleimputation("H7", queryalignments)
#exit()

#segments = AbstractString["H5HA","H5PB1", "H5PB2", "H5M", "H5NS", "H5PA", "H5NA", "H5NP"] 
#queryalignments = AbstractString["../data/Metric_Datasets/H5/HA/H5_HA_sampled_P.fasta", "../data/Metric_Datasets/H5/PB1/H5_PB1_sampled_P.fasta", "../data/Metric_Datasets/H5/PB2/H5_PB2_sampled_P.fasta", "../data/Metric_Datasets/H5/M/H5_M_sampled_P.fasta", "", ""]
#queryalignments = AbstractString["../data/Metric_Datasets/H5/HA/H5_HA_P.fasta", "../data/Metric_Datasets/H5/PB1/H5_PB1_P.fasta", "../data/Metric_Datasets/H5/PB2/H5_PB2_P.fasta", "../data/Metric_Datasets/H5/M/H5_M_P.fasta", "", "../data/Metric_Datasets/H5/PA/PA_short_P.fasta", "../data/Metric_Datasets/H5/PA/PA_short_P.fasta", "../data/Metric_Datasets/H5/NA/NA_H5NX_short_F_P.fasta", "../data/Metric_Datasets/H5/NP/NP_short_P.fasta"]
#ParallelEvolution.creatematrix("H5",queryalignments)
#ParallelEvolution.creatematrix("H5")
#ParallelEvolution.simpleimputation("H5",queryalignments)
#ParallelEvolution.simpleimputation("H5")
#exit()
#mappingsequence = "MNTQILAFIACMLIGTKGDKICLGHHAVANGTKVNTLTERGIEVVNATETVETVNIKKICTQGKRPTDLGQCGLLGTLIGPPQCDQFLEFDANLIIERREGTDVCYPGKFTNEESLRQILRGSGGIDKESMGFTYSGIRTNGATSACRRSGSSFYAEMKWLLSNSDNAAFPQMTKSYRNPRNKPALIIWGVHHSGSATEQTKLYGSGNKLITVGSSKYQQSFTPSPGARPQVNGQSGRIDFHWLLLDPNDTVTFTFNGAFIAPDRASFFRGESLGVQSDVPLDSGCEGDCFHSGGTIVSSLPFQNINPRTVGKCPRYVKQTSLLLATGMRNVPENPKRGLFGAIAGFIENGWEGLIDGWYGFRHQNAQGEGTAADYKSTQSAIDQITGKLNRLIDKTNQQFELIDNEFSEIEQQIGNVINWTRDSMTEVWSYNAELLVAMENQH-------------------------------------------------------------------------------------------------------------------"
#ParallelEvolution.scoresequence(mappingsequence)
#println(ParallelEvolution.marginaltraitprobability(1, 1, 0.01, 0.9))
#using LinearAlgebra
#Q,freqs = ParallelEvolution.getIndelTmatrix(0.5, 0.25, 1.0, 0.4, 10.0)
#=
println("DEL HP ", freqs[(1-1)*2 + 1])
println("DEL LP ", freqs[(1-1)*2 + 2])
println("INS HP ", freqs[(2-1)*2 + 1])
println("INS LP ", freqs[(2-1)*2 + 2])
println(Q)
println(exp(Q*100.0))
println(freqs)
println(transpose(freqs)*Q)
#println(ParallelEvolution.getLGmatrixwithdeletion(0.2, 2.0))
=#
#exit()

parsed_args = ParallelEvolution.parse_commandline()
if parsed_args["mlmodel"]
    ParallelEvolution.mlmodel()
elseif parsed_args["mlindelmodel"]
    ParallelEvolution.mlindelmodel()
elseif parsed_args["score"]
    ParallelEvolution.score()
elseif parsed_args["coevolution"]
    ParallelEvolution.coevolution()
elseif parsed_args["simulation"]
    ParallelEvolution.simulation()
elseif parsed_args["processsimulations"]
    ParallelEvolution.processsimulations()
else
    ParallelEvolution.bhattacharya()
end
