module LG

function load_parameters(parameter_file)
  f = open(parameter_file)
  lines = readlines(f)
  close(f)

  S = zeros(Float64,20, 20)


  for i=1:20
   spl = split(lines[i])
   for j=1:length(spl)
     S[i+1,j] = parse(Float64, spl[j])
     S[j,i+1] = S[i+1,j]
   end
  end

  eqfreqs = zeros(Float64,20)
  spl = split(lines[21])
  for i=1:20
    eqfreqs[i] = parse(Float64, spl[i])
  end

  return S,eqfreqs
end

export LGmatrix,LGfreqs
S, LGfreqs = load_parameters(joinpath(@__DIR__, "lg_LG.PAML.txt"))
LGmatrix = zeros(20,20)
for i=1:20
    for j=1:20
        if i != j
            LGmatrix[i,j] = S[i,j]*LGfreqs[j]
            LGmatrix[i,i] -= LGmatrix[i,j]
        end
    end
end

end
