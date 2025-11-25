"""
`TimeFilter(callbacks...; times)`

A filter that calls each filter in `callbacks` once per each time in `times`. 
This is intended to run in the innest cycle of `run_RD` and it's very fast
"""
function TimeFilter(callbacks...; times)
    next = Ref(1)
    function stats(t, s)
        while next[] <= lastindex(times) && times[next[]] <= t
            for cb in callbacks
                cb(times[next[]], s)
            end
            next[] += 1
        end
    end
end

"""
`StopWatchFilter(callbacks...; seconds=1)`

A filter that calls `callbacks` with at least `seconds` wall-clock 
delay between calls (intended e.g. for iterative plotting)
"""
function StopWatchFilter(callbacks...; seconds=1)
    last = time()
    function stats(t, s)
        next = time()
        if next - last > seconds
            last = next
            for cb in callbacks
                cb(t, s)
            end
        end
    end
end
