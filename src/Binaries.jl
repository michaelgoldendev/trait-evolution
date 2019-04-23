module Binaries
    using FastaIO
    using SHA
    push!(LOAD_PATH,@__DIR__)
    using CommonUtils

    figtree_binary  = joinpath(@__DIR__,"..","binaries","figtree-custom.jar")
    function getfigtreesvg(newick_file::AbstractString, width::Int=320, height::Int=500)
        return read(`java -jar $(figtree_binary) -graphic SVG -width $(width) -height $(height) $(newick_file)`, String)
    end

    function fasttreegtr(alignmentfile::AbstractString)
        fastastring = open(alignmentfile) do file
            read(file, String)
        end
        cachepath = joinpath(@__DIR__,"..","cache")
        mkpath(cachepath)
        cachefile = joinpath(cachepath, string(CommonUtils.sha256base36(fastastring), ".fasttreegtr.nwk"))
        if isfile(cachefile)
            newickstring = open(cachefile) do file
                read(file, String)
            end
            newickstring, cachefile
        else
            fasttreepath = fasttreepath = joinpath(@__DIR__,"..","binaries","FastTreeMP")
            if Sys.iswindows()
              fasttreepath = joinpath(@__DIR__,"..","binaries","FastTree.exe")
            end
            newickstring = read(`$fasttreepath -nt -gtr -nosupport $alignmentfile`, String)

            fout = open(cachefile, "w")
            print(fout, strip(newickstring))
            close(fout)
            return newickstring, cachefile
        end
    end

    function fasttreeaa(alignmentfile::AbstractString)
        fastastring = open(alignmentfile) do file
            read(file, String)
        end
        cachepath = joinpath(@__DIR__,"..","cache")
        mkpath(cachepath)
        cachefile = joinpath(cachepath, string(CommonUtils.sha256base36(fastastring), ".fasttreeaa.nwk"))
        if isfile(cachefile)
            newickstring = open(cachefile) do file
                read(file, String)
            end
            newickstring, cachefile
        else
            fasttreepath = fasttreepath = joinpath(@__DIR__,"..","binaries","FastTreeMP")
            if Sys.iswindows()
              fasttreepath = joinpath(@__DIR__,"..","binaries","FastTree.exe")
            end
            newickstring = read(`$fasttreepath -nosupport $alignmentfile`, String)

            fout = open(cachefile, "w")
            print(fout, strip(newickstring))
            close(fout)
            return newickstring, cachefile
        end
    end
end
