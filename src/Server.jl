using Distributed
addprocs(2)
@everywhere push!(LOAD_PATH,@__DIR__)
@everywhere using TaskScheduling
using Binaries
using CommonUtils
using Mux
using Mustache
using HttpCommon



function serveassets(subpath, req)
    filepath = joinpath(subpath, req[:params][:file])
    ret = ""
    if isfile(filepath)
        ret = open(filepath) do file
            read(file, String)
        end
        if endswith(filepath, ".css")
            headers  = HttpCommon.headers()
            headers["Content-Type"] = "text/css"
            return Dict(
            :headers => headers,
            :body=> ret
            )
        elseif endswith(filepath, ".js")
            headers  = HttpCommon.headers()
            headers["Content-Type"] = "text/js"
            return Dict(
            :headers => headers,
            :body=> ret
            )
        elseif endswith(filepath, ".html")
            headers  = HttpCommon.headers()
            headers["Content-Type"] = "text/html"
            return Dict(
            :headers => headers,
            :body=> ret
            )
        elseif endswith(filepath, ".png")
            headers  = HttpCommon.headers()
            headers["Content-Type"] = "image/png"
            return Dict(
            :headers => headers,
            :body=> ret
            )
        else
            return ret
        end
    end

    return HttpCommon.Response(404)
end

function redirect(req, scheduler)
    println("000")
    queryparams = HttpCommon.parsequerystring(req[:query])
    taskid = queryparams["taskid"]
    println("AAA")
    taskstatus = gettaskstatus(scheduler, taskid)
    println("BBB")
    value = fetchvalue(scheduler, taskid)
    println("CCC")

    headers  = HttpCommon.headers()
    headers["Content-Type"] = "application/json"
    if value != nothing
        println("DDD")
        body = "{\"taskid\": $(taskid), \"taskstatus\": $(taskstatus), \"downloadurl\": \"$(value)\"}"
        println(body)
        return Dict(
            :headers => headers,
            :body=> body
        )
    else
        println("EEE")
        return Dict(
            :headers => headers,
            :body=> "{\"taskid\": $(taskid), \"taskstatus\": $(taskstatus)}"
        )
    end
end

function uploadevent(req, scheduler)
    datastring = read(Base.IOBuffer(req[:data]),String)
    filepath = joinpath(@__DIR__,"..","cache", string(CommonUtils.sha256base36(datastring),".alignment.fas"))
    fout = open(filepath, "w")
    write(fout, datastring)
    close(fout)

    outfile =  string("cache/", CommonUtils.sha256base36(datastring),".bhattacharya.csv")
    t1() = TaskScheduling.runbhattacharya(filepath, outfile)
    task = CustomTask(t1)
    taskid = addtask(scheduler, task)

    return Dict(
    :headers => HttpCommon.headers(),
    :body=> "{\"taskid\": $(taskid)}"
    )
end

scheduler = TaskScheduler(2)

@app test = (
  Mux.defaults,
  route("/upload/alignment", req -> uploadevent(req, scheduler)),
  page("/redirect", req -> redirect(req, scheduler)),
  page("/:file", req -> serveassets("assets/",req)),
  page("/cache/:file", req -> serveassets("cache/",req)),
  page("/css/:file", req -> serveassets("assets/css",req)),
  page("/data/:file", req -> serveassets("assets/data",req)),
  page("/fonts/:file", req -> serveassets("assets/fonts",req)),
  page("/images/:file", req -> serveassets("assets/images",req)),
  page("/js/:file", req -> serveassets("assets/js",req)),
  Mux.notfound())

serve(test)
loop(scheduler)
#=
println("Waiting for 'exit()' command.\n")
for line in eachline(stdin)
    cmdstr = strip(line)
    if cmdstr == "exit()"
      println("Exiting...")
      break
    else
      println("Command '$cmdstr' not recognized.\n")
    end
end
exit()
=#
