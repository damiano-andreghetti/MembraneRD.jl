using ColorVectorSpace, Images, Colors, Compose
import Cairo, Fontconfig

function hexagon(x, y, r)
    polygon([[(x[i]+r*cos(π/3*(k+3/2)),
               y[i]+r*sin(π/3*(k+3/2))) for k in 0:6] 
                    for i in eachindex(x,y)])
end

"""
`composed(s::State; posx, posy, colors, Δ)`

Generates a `Compose` image of an hexagonal grid from `s::State` 
"""
function composed(s::State; posx, posy, colors, Δ)
    X, Y = Δ .* posx, Δ .* posy
    (x0,x1),(y0,y1) = extrema(X), extrema(Y)
    set_default_graphic_size(x1-x0+3Δ, y1-y0+3Δ)
    compose(context(), 
        hexagon(X, Y, 1Δ),
        fill(s.membrane * colors)
    ) 
end

"""
`Plotter(posx, posy; colors, Δ=mm)`

A filter to generate a `Compose` image from each state `s`. Use e.g. the
`display ∘ Plotter(posx, posy; colors)` filter for iterative display, etc.
"""
function Plotter(posx, posy; colors, Δ=mm)
    f(_, s::State) = composed(s; posx, posy, colors, Δ)
end

"""
`raster(c::Compose.Context) -> Matrix{<:RGB}`

Converts a `Context` object into a raster image.
"""
function raster(c)
    io = IOBuffer()
    draw(PNG(io; emit_on_finish=false), c)
    A = Images.load(io)
    m, n = 2 .* (size(A) .÷ 2)
    #crop to even dimensions and N0f8 colorspace
    RGB{N0f8}.(@view A[1:m,1:n])
end
