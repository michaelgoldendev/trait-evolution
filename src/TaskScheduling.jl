module TaskScheduling
    using Distributed

    push!(LOAD_PATH,@__DIR__)
    using ParallelEvolution

    function runbhattacharya(fastafile, outfile)
        try
            println(fastafile)
            newickstring, treefile = Binaries.fasttreeaa(fastafile)
            return ParallelEvolution.bhattacharya(fastafile, treefile, split(".LP.,.HP.",","), outfile)
        catch e
            return nothing
        end
        return nothing
    end

    export CustomTask
    mutable struct CustomTask
        #lock::ReentrantLock
        func::Function
        running::Bool
        finished::Bool
        future::Future
        ret::Any

        function CustomTask(func)
            new(func,false, false, Future(),nothing)
        end
    end

    export istaskdone
    function istaskdone(task::CustomTask)
        finished = false
        #Base.lock(task.lock)
        finished = Distributed.isready(task.future)
        #Base.unlock(task.lock)
        return finished
    end

    export schedule
    function schedule(task::CustomTask)
        #Base.lock(task.lock)
        task.running = true
        task.future = Distributed.@spawn task.func()
        #Base.unlock(task.lock)
        #=
        @async begin

        end=#
    end

    export fetch
    function fetch(task::CustomTask)
        ret = Distributed.fetch(task.future)
        return ret
    end


    export TaskScheduler
    mutable struct TaskScheduler
        lock::ReentrantLock
        queuedtasks::Dict{AbstractString,CustomTask}
        activetasks::Dict{AbstractString,CustomTask}
        completedtasks::Dict{AbstractString,CustomTask}
        maxactivetasks::Int
        taskcount::Int

        function TaskScheduler(maxactivetasks::Int)
            new(ReentrantLock(), Dict{AbstractString,CustomTask}(), Dict{AbstractString,CustomTask}(), Dict{AbstractString,CustomTask}(), maxactivetasks, 0)
        end
    end

    export addtask
    function addtask(scheduler::TaskScheduler, task::CustomTask)
        Base.lock(scheduler.lock)
        taskid = string(scheduler.taskcount)
        scheduler.queuedtasks[taskid] = task
        scheduler.taskcount += 1
        Base.unlock(scheduler.lock)
        update(scheduler)
        return taskid
    end

    export gettaskstatus
    function gettaskstatus(scheduler::TaskScheduler, taskid::AbstractString)
        update(scheduler)
        ret = 0
        Base.lock(scheduler.lock)
        if haskey(scheduler.completedtasks, taskid)
            ret = 3
        elseif haskey(scheduler.activetasks, taskid)
            ret = 2
        elseif haskey(scheduler.queuedtasks, taskid)
            ret = 1
        end
        Base.unlock(scheduler.lock)
        return ret
    end



    export fetchvalue
    function fetchvalue(scheduler::TaskScheduler, taskid::AbstractString)
        update(scheduler)
        println("fetching $(taskid)")
        ret = nothing
        Base.lock(scheduler.lock)
        if haskey(scheduler.completedtasks, taskid) && istaskdone(scheduler.completedtasks[taskid])
            println("key ", [k for k in keys(scheduler.completedtasks)])
            ret = fetch(scheduler.completedtasks[taskid])
        end
        Base.unlock(scheduler.lock)
        println("finished fetching $(taskid)")
        return ret
    end

    function update(scheduler::TaskScheduler)
        println("W8")
        println("W0")
        println("queued: $(length(scheduler.queuedtasks)), active: $(length(scheduler.activetasks))")
        Base.lock(scheduler.lock)

        println("W1")
        keylist = [k for k in keys(scheduler.activetasks)]
        for taskid in keylist
            if istaskdone(scheduler.activetasks[taskid])
                scheduler.completedtasks[taskid] = scheduler.activetasks[taskid]
                delete!(scheduler.activetasks, taskid)
            end
        end
        println("W2")

        keylist = [k for k in keys(scheduler.queuedtasks)]
        for taskid in keylist
            println("W3")
            if length(scheduler.activetasks) < scheduler.maxactivetasks
                task = scheduler.queuedtasks[taskid]
                delete!(scheduler.queuedtasks, taskid)
                scheduler.activetasks[taskid] = task
                println("W4")
                schedule(task)
                println("W5")
            end
        end
        println("W6")
        Base.unlock(scheduler.lock)
        println("W7")
    end

    export loop
    function loop(scheduler::TaskScheduler)
        while true
            update(scheduler)
            sleep(10.0)
        end
    end
end
