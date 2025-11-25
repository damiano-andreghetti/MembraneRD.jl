import VideoIO

"""
`savevideo(fn, stack, plotter; kwd)`

Saves a .MP4 file generated from a stack `State`s 

`fn`: filename
`stack::AbstractVector`: stack of images
"""
function savevideo(fn, stack; framerate=10, encoder_options = (crf=23, preset="medium"), showprogress=true)
    VideoIO.open_video_out(fn, first(stack); framerate, encoder_options) do w
        if showprogress
            @showprogress for img in stack
                write(w, img)
            end
        else
            for img in stack
                write(w, img)
            end
        end
    end
#    VideoIO.save(fn, stack; framerate, encoder_options)
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
    VideoIO.open_video_out(w->f((t,s)->write(w,plotter(t,s))), fn, plotter(0.0,s); framerate, encoder_options)
end