"""
Pushes something into a vector
"""
function Pusher(::Type{T} = Any) where T
    stack = T[]
    pusher(x...) = push!(stack, deepcopy(x))
end

