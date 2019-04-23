module ChiSquaredEvolution
    using MolecularEvolution
    using CommonUtils
    using HypothesisTests
    using Distributions

    export expectedchi
    function expectedchi(observed::Array{Float64,2})
        chi2 = 0.0
        expected = zeros(Float64, size(observed,1), size(observed,2))
        total = sum(observed)
        for i=1:size(observed,1)
            for j=1:size(observed,2)
                expected[i,j] = sum(observed[i,:])*sum(observed[:,j])/total
                chi2 += ((observed[i,j] - expected[i,j])^2.0)/expected[i,j]
            end
        end
        return expected,chi2
    end

    export fisherstest
    function fisherstest(observed::Array{Float64,2})
        #=
        chi2 = 0.0
        expected = zeros(Float64, size(observed,1), size(observed,2))
        total = sum(observed)
        for i=1:size(observed,1)
            for j=1:size(observed,2)
                expected[i,j] = sum(observed[i,:])*sum(observed[:,j])/total
                chi2 += ((observed[i,j] - expected[i,j])^2.0)/expected[i,j]
            end
        end

        D = abs(chi2)
        if D <= 0.0
                return 1.0
        end
        return ccdf(Chisq(1), D)=#

        pval = 1.0
        try
            pval = pvalue(FisherExactTest(convert(Int,observed[1,1]),convert(Int,observed[1,2]),convert(Int,observed[2,1]),convert(Int,observed[2,2])))
        catch

        end
        return pval
    end

    export chisquaredtest
    function chisquaredtest(nodelist::Array{TreeNode,1}, markednodes::Array{Int,1}, traits::Array{Int,1}, samples::Array{Int,2}, targetaa::Int, collapseall::Bool=false)
        nsamples = size(samples,1)
        chi2_res1 = Float64[]
        chi2_res2 = Float64[]
        chi2_res3 = Float64[]
        fisher_res1 = Float64[]
        fisher_res2 = Float64[]
        fisher_res3 = Float64[]
        for s=1:nsamples
            observed1 = zeros(Float64,2,2)
            observed2 = zeros(Float64,2,2)
            observed3 = zeros(Float64,2,2)
            for nodeindex=1:length(markednodes)
                node = nodelist[nodeindex]
                if (collapseall && markednodes[nodeindex] > 0) || ((!collapseall && markednodes[nodeindex] == 2) || (!collapseall && isleafnode(node) && traits[node.seqindex] == 1))
                    parentaa = ((samples[s,get(nodelist[nodeindex].parent).nodeindex]-1) % 20) + 1
                    nodeaa = ((samples[s,nodeindex]-1) % 20) + 1
                    trait = div(samples[s,nodeindex]-1,20) + 1
                    if (collapseall && markednodes[nodeindex] > 0)
                        trait = markednodes[nodeindex]
                    elseif !collapseall && markednodes[nodeindex] == 2
                        trait = markednodes[nodeindex]
                    elseif !collapseall && isleafnode(node) && traits[node.seqindex] == 1
                        trait = traits[node.seqindex]
                    end

                    if parentaa == targetaa && nodeaa != targetaa
                        observed1[1,trait] += 1.0
                    elseif parentaa == targetaa && nodeaa == targetaa
                        observed1[2,trait] += 1.0
                    end

                    if parentaa != targetaa && nodeaa == targetaa
                        observed2[2,trait] += 1.0
                    elseif parentaa != targetaa && nodeaa != targetaa
                        observed2[1,trait] += 1.0
                    end

                    if nodeaa == targetaa
                        observed3[2,trait] += 1.0
                    elseif nodeaa != targetaa
                        observed3[1,trait] += 1.0
                    end
                end
            end

            #println(observed1)
            expected1, chi2_1 = expectedchi(observed1)
            #println(expected1)
            #println(chi2_1)
            chi2_1reported = -chi2_1
            if observed1[2,2] > expected1[2,2]
                chi2_1reported = chi2_1
            end
            push!(chi2_res1, chi2_1reported)
            push!(fisher_res1, fisherstest(observed1))
            #println(chi2_1reported)
            #println()

            #println(observed2)
            expected2, chi2_2 = expectedchi(observed2)
            #println(expected2)
            #println(chi2_2)
            chi2_2reported = -chi2_2
            if observed2[2,2] > expected2[2,2]
                chi2_2reported = chi2_2
            end
            push!(chi2_res2, chi2_2reported)
            push!(fisher_res2, fisherstest(observed2))
            #println(chi2_2reported)
            #println()

            #println(observed3)
            expected3, chi2_3 = expectedchi(observed3)
            #println(expected3)
            #println(chi2_3)
            chi2_3reported = -chi2_3
            if observed3[2,2] > expected3[2,2]
                chi2_3reported = chi2_3
            end
            push!(chi2_res3, chi2_3reported)
            push!(fisher_res3, fisherstest(observed3))
            #println(chi2_3reported)
            #println()
        end
        return safemedian(chi2_res1),safepercentile(chi2_res1,5.0),safepercentile(chi2_res1,95.0), safemedian(chi2_res2), safepercentile(chi2_res2,5.0), safepercentile(chi2_res2,95.0), safemedian(chi2_res3), safepercentile(chi2_res3,5.0), safepercentile(chi2_res3,95.0),safemedian(fisher_res1,1.0),safemedian(fisher_res2,1.0),safemedian(fisher_res3,1.0)
    end

    export chisquaredtestmedians
    function chisquaredtestmedians(nodelist::Array{TreeNode,1}, markednodes::Array{Int,1}, traits::Array{Int,1}, samples::Array{Int,2}, targetaa::Int, collapseall::Bool=false)
        nsamples = size(samples,1)
        chi2_res1 = Float64[]
        chi2_res2 = Float64[]
        chi2_res3 = Float64[]
        fisher_res1 = Float64[]
        fisher_res2 = Float64[]
        fisher_res3 = Float64[]
        for s=1:nsamples
            observed1 = zeros(Float64,2,2)
            observed2 = zeros(Float64,2,2)
            observed3 = zeros(Float64,2,2)
            for nodeindex=1:length(markednodes)
                node = nodelist[nodeindex]
                if (collapseall && markednodes[nodeindex] > 0) || ((!collapseall && markednodes[nodeindex] == 2) || (!collapseall && isleafnode(node) && traits[node.seqindex] == 1))
                    parentaa = ((samples[s,get(nodelist[nodeindex].parent).nodeindex]-1) % 20) + 1
                    nodeaa = ((samples[s,nodeindex]-1) % 20) + 1
                    trait = div(samples[s,nodeindex]-1,20) + 1
                    if (collapseall && markednodes[nodeindex] > 0)
                        trait = markednodes[nodeindex]
                    elseif !collapseall && markednodes[nodeindex] == 2
                        trait = markednodes[nodeindex]
                    elseif !collapseall && isleafnode(node) && traits[node.seqindex] == 1
                        trait = traits[node.seqindex]
                    end

                    if parentaa == targetaa && nodeaa != targetaa
                        observed1[1,trait] += 1.0
                    elseif parentaa == targetaa && nodeaa == targetaa
                        observed1[2,trait] += 1.0
                    end

                    if parentaa != targetaa && nodeaa == targetaa
                        observed2[2,trait] += 1.0
                    elseif parentaa != targetaa && nodeaa != targetaa
                        observed2[1,trait] += 1.0
                    end

                    if nodeaa == targetaa
                        observed3[2,trait] += 1.0
                    elseif nodeaa != targetaa
                        observed3[1,trait] += 1.0
                    end
                end
            end

            #println(observed1)
            expected1, chi2_1 = expectedchi(observed1)
            #println(expected1)
            #println(chi2_1)
            chi2_1reported = -chi2_1
            if observed1[2,2] > expected1[2,2]
                chi2_1reported = chi2_1
            end
            push!(chi2_res1, chi2_1reported)
            push!(fisher_res1, fisherstest(observed1))
            #println(chi2_1reported)
            #println()

            #println(observed2)
            expected2, chi2_2 = expectedchi(observed2)
            #println(expected2)
            #println(chi2_2)
            chi2_2reported = -chi2_2
            if observed2[2,2] > expected2[2,2]
                chi2_2reported = chi2_2
            end
            push!(chi2_res2, chi2_2reported)
            push!(fisher_res2, fisherstest(observed2))
            #println(chi2_2reported)
            #println()

            #println(observed3)
            expected3, chi2_3 = expectedchi(observed3)
            #println(expected3)
            #println(chi2_3)
            chi2_3reported = -chi2_3
            if observed3[2,2] > expected3[2,2]
                chi2_3reported = chi2_3
            end
            push!(chi2_res3, chi2_3reported)
            push!(fisher_res3, fisherstest(observed3))
            #println(chi2_3reported)
            #println()
        end
        return safemedian(chi2_res1), safemedian(chi2_res2), safemedian(chi2_res3), safemedian(fisher_res1,1.0),safemedian(fisher_res2,1.0),safemedian(fisher_res3,1.0)
    end

    export chisquaredtestmedians
    function chisquaredtestmedians(nodelist::Array{TreeNode,1}, markednodes::Array{Int,1}, samples::Array{Int,2}, targetaa::Int)
        nsamples = size(samples,1)
        chi2_res1 = Float64[]
        chi2_res2 = Float64[]
        chi2_res3 = Float64[]
        fisher_res1 = Float64[]
        fisher_res2 = Float64[]
        fisher_res3 = Float64[]
        for s=1:nsamples
            observed1 = zeros(Float64,2,2)
            observed2 = zeros(Float64,2,2)
            observed3 = zeros(Float64,2,2)
            for nodeindex=1:length(markednodes)
                node = nodelist[nodeindex]
                trait = markednodes[nodeindex]
                if trait > 0
                    parentaa = ((samples[s,get(nodelist[nodeindex].parent).nodeindex]-1) % 20) + 1
                    nodeaa = ((samples[s,nodeindex]-1) % 20) + 1
                    if parentaa == targetaa && nodeaa != targetaa
                        observed1[1,trait] += 1.0
                    elseif parentaa == targetaa && nodeaa == targetaa
                        observed1[2,trait] += 1.0
                    end

                    if parentaa != targetaa && nodeaa == targetaa
                        observed2[2,trait] += 1.0
                    elseif parentaa != targetaa && nodeaa != targetaa
                        observed2[1,trait] += 1.0
                    end

                    if nodeaa == targetaa
                        observed3[2,trait] += 1.0
                    elseif nodeaa != targetaa
                        observed3[1,trait] += 1.0
                    end
                end
            end

            #println(observed1)
            expected1, chi2_1 = expectedchi(observed1)
            #println(expected1)
            #println(chi2_1)
            chi2_1reported = -chi2_1
            if observed1[2,2] > expected1[2,2]
                chi2_1reported = chi2_1
            end
            push!(chi2_res1, chi2_1reported)
            push!(fisher_res1, fisherstest(observed1))
            #println(chi2_1reported)
            #println()

            #println(observed2)
            expected2, chi2_2 = expectedchi(observed2)
            #println(expected2)
            #println(chi2_2)
            chi2_2reported = -chi2_2
            if observed2[2,2] > expected2[2,2]
                chi2_2reported = chi2_2
            end
            push!(chi2_res2, chi2_2reported)
            push!(fisher_res2, fisherstest(observed2))
            #println(chi2_2reported)
            #println()

            #println(observed3)
            expected3, chi2_3 = expectedchi(observed3)
            #println(expected3)
            #println(chi2_3)
            chi2_3reported = -chi2_3
            if observed3[2,2] > expected3[2,2]
                chi2_3reported = chi2_3
            end
            push!(chi2_res3, chi2_3reported)
            push!(fisher_res3, fisherstest(observed3))
            #println(chi2_3reported)
            #println()
        end
        return safemedian(chi2_res1), safemedian(chi2_res2), safemedian(chi2_res3), safemedian(fisher_res1,1.0),safemedian(fisher_res2,1.0),safemedian(fisher_res3,1.0)
    end
end
