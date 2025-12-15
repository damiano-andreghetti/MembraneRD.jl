"""
Pushes something into a vector
"""
function Pusher(f, ::Type{T} = Any) where T
    stack = T[]
    pusher(xs...) = push!(stack, f(xs...))
end

