using ProgressMeter

"""
`ProgressShower(T; kwds...)`

A filter to show progress

`T`: final time
"""
function ProgressShower(T; kwds...)
    p = Progress(100; kwds...)
    stat(t, _) = update!(p, floor(Int, 100*t/T))
end