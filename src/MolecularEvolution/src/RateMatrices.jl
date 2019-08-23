export gtr
function gtr(q1::Float64, q2::Float64, q3::Float64, q4::Float64, q5::Float64, q6::Float64, piv::Array{Float64,1}=Float64[0.25,0.25,0.25,0.25])
    S = zeros(Float64,4,4)
    S[1,2] = q1
    S[1,3] = q2
    S[1,4] = q3

    S[2,1] = q1
    S[2,3] = q4
    S[2,4] = q5

    S[3,1] = q2
    S[3,2] = q4
    S[3,4] = q6

    S[4,1] = q3
    S[4,2] = q5
    S[4,3] = q6

    for i=1:size(S,1)
        for j=1:size(S,2)
            S[i,j] *= piv[j]
        end
    end

    Q = copy(S)
    for i=1:size(Q,1)
        diagelem = 0.0
        for j=1:size(Q,2)
            if i != j
                diagelem -= S[i,j]
            end
        end
        Q[i,i] = diagelem
    end
    return Q
end

export calculatemusefreqs
function calculatemusefreqs(obsFreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64)
    dinucfreqs = zeros(Float64,16)

    piGC = obsFreqs[2]*obsFreqs[3]
    piAT = obsFreqs[1]*obsFreqs[4]
    piGT = obsFreqs[3]*obsFreqs[4]

    kappa = 1.0/(1.0 + 2.0*(piAT*(lambdaAT*lambdaAT-1.0)) + 2.0*(piGC*(lambdaGC*lambdaGC-1.0)) + 2.0*(piGT*(lambdaGT*lambdaGT-1.0)))

    basepairings = Int[0,0,0,2,0,0,1,0,0,1,0,3,2,0,3,0]

    for h=1:4
        for v=1:4
            idx = (h-1)*4+v
            if basepairings[idx] == 1
                dinucfreqs[idx] = kappa*lambdaGC*lambdaGC*obsFreqs[h]*obsFreqs[v]
            elseif basepairings[idx] == 2
                dinucfreqs[idx] = kappa*lambdaAT*lambdaAT*obsFreqs[h]*obsFreqs[v]
            elseif basepairings[idx] == 3
                dinucfreqs[idx] = kappa*lambdaGT*lambdaGT*obsFreqs[h]*obsFreqs[v]
            else
                dinucfreqs[idx] = kappa*obsFreqs[h]*obsFreqs[v]
            end
        end
    end
    return dinucfreqs
end

export muse
function muse(lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64,q1::Float64, q2::Float64, q3::Float64, q4::Float64, q5::Float64, q6::Float64, obsfreqs::Array{Float64,1}=Float64[0.25,0.25,0.25,0.25], siteRate1::Float64=1.0, siteRate2::Float64=1.0)
    gtrQ = gtr(q1,q2,q3,q4,q5,q6,ones(Float64,4))
    basepairings = Int[0,0,0,2,0,0,1,0,0,1,0,3,2,0,3,0]
    dinucfreqs = calculatemusefreqs(obsfreqs,lambdaGC,lambdaAT,lambdaGT)
    Q = zeros(Float64, 16, 16)
    for h=1:16
        for v=h+1:16
            if v != h
                fromNuc = 0
                toNuc   = 0
                if div(h-1,4)+1 == div(v-1,4)+1
                    toNuc   = ((v-1) % 4) + 1
                    fromNuc = ((h-1) % 4) + 1
                elseif ((v-1)%4) + 1 == ((h-1) %4)+1
                    toNuc   = div(v-1,4) + 1
                    fromNuc = div(h-1,4) + 1
                end
                if fromNuc > 0
                    rateMult  = 1.0
                    rateMult2 = 1.0

                    if basepairings[v] == 1
                        rateMult  *= lambdaGC
                        rateMult2 /= lambdaGC
                    end
                    if basepairings[h] == 1
                        rateMult /= lambdaGC
                        rateMult2 *= lambdaGC
                    end


                    if basepairings[v] == 2
                        rateMult  *= lambdaAT
                        rateMult2 /= lambdaAT
                    end
                    if basepairings[h] == 2
                        rateMult /= lambdaAT
                        rateMult2 *= lambdaAT
                    end

                    if basepairings[v] == 3
                        rateMult  *= lambdaGT
                        rateMult2 /= lambdaGT
                    end
                    if basepairings[h] == 3
                        rateMult /= lambdaGT
                        rateMult2 *= lambdaGT
                    end

                    if fromNuc != toNuc
                        rateMult  *= gtrQ[fromNuc,toNuc]
                        rateMult2 *= gtrQ[fromNuc,toNuc]
                    end

                    if div(h-1,4) == div(v-1,4)
                        rateMult  *= siteRate2
                        rateMult2 *= siteRate2
                    elseif (v-1) %4  == (h-1) %4
                        rateMult  *= siteRate1
                        rateMult2 *= siteRate1
                    end

                    Q[h,v] = rateMult*obsfreqs[toNuc]
                    Q[v,h] = rateMult2*obsfreqs[fromNuc]
                end
            end
        end
    end

    for i=1:16
        Q[i,i] = 0.0
        for j=1:16
            if i != j
                Q[i,i] -= Q[i,j]
            end
        end
    end

    return Q
end
