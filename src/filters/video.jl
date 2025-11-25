import VideoIO

"""
`savevideo(fn, stack, plotter; kwd)`

Saves a .MP4 file generated from a stack of either `State`s or `(::Float64,::State)` tuples 

`fn`: filename
`stack::AbstractVector`: stack of images
`plotter`: converts a `t,s` into an image (e.g. a `Plotter` works)
"""
function savevideo(fn, stack::AbstractVector{Tuple{Float64,State}}, plotter; framerate=10, encoder_options = (crf=23, preset="medium"))
    savevideo(fn, [plotter(t, s) for (t,s) in stack]; framerate, encoder_options)
end

function savevideo(fn, stack::AbstractVector{State}, plotter; framerate=10, encoder_options = (crf=23, preset="medium"))
    savevideo(fn, [plotter(0.0, s) for s in stack]; framerate, encoder_options)
end

"""
`savevideo(fn, stack, plotter; kwd)`

Saves a .MP4 file generated from a stack of `RGB` images or `Compose.Context` images

`fn`: filename
`stack::AbstractVector`: stack of images
"""
function savevideo(fn, stack::AbstractVector{<:Compose.Context}; framerate=10, encoder_options = (crf=23, preset="medium"))
    savevideo(fn, raster.(stack); framerate, encoder_options)
end

function savevideo(fn, stack::AbstractVector{<:AbstractMatrix}; framerate=10, encoder_options = (crf=23, preset="medium"))
    VideoIO.save(fn, stack; framerate, encoder_options)
end

"""
`OpenVideo(f, filename, s::State; plotter, kwd...)`

This function is intended for iterative video generation, without saving
intermediate raster images and thus saving memory/disk space

use e.g. as follows:

```julia
OpenVideo("video.mp4", s; plotter) do framesaver
    stats = TimeFilter(framesaver; times=Tmeas:Tmeas:T)
    run_RD!(s, M, T; stats, rng)
end
```

where `plotter` takes `t,s` and produces a bitmap matrix.
"""
function OpenVideo(f, fn, s::State; plotter, framerate=10, encoder_options = (crf=23, preset="medium"))
    VideoIO.open_video_out(w->f((t,s)->write(w,plotter(t,s))), fn, plotter(0.,s); framerate, encoder_options)
end