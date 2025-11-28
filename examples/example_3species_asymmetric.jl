using MembraneRD, Random, Colors, MembraneRD.Filters

species = @species A B C EBA ECA EAB EAC

function build_model3(L)
    G, layout = MembraneRD.gen_hex_lattice(L)
    N = nv(G)
    cat = (@catalytic (1.0,1.0) EBA + B => EBA + A),
        (@catalytic (1.0,1.0) EAB + A => EAB + B),
        (@catalytic (1.0,1.0) EAC + A => EAC + C),
        (@catalytic (1.0,1.0) ECA + C => ECA + A)
    att = (EBA,A,1.0), (EAB,B,1.0), (ECA,A,1.0), (EAC,C,1.0)
    det = (EBA,1.0), (EAB,1.0), (ECA,1.0), (EAC,1.0)
    dif = (A,0.1), (B,0.1), (C,0.1), 
        (EBA,0.1), (ECA,0.1), (EAB,0.1), (EAC,0.1)
    M = Model(; species, G, cat, att, det, dif, rho_0 = 0.0)
    (; M, mem = [2N,10N,10N,0,0,0,0], 
        cyto = floor.(Int,[0,0,0,0.3N,0.3N,0.5N,0.5N]), 
        layout)
end

T = 20000.0
L = 150

#for reproducibility
rng = Random.Xoshiro(22)
(; M, mem, cyto, layout) = build_model3(L)
s = State(M, mem, cyto; rng)

times = 0:T/200:T
saver = Pusher((t,s,_...)->(t,deepcopy(s)), Tuple{Float64,State})
colors = [RGB(m == A, m == B, m == C) for m in species]/30
plotter = Plotter(layout; colors)
displayer = StopWatchFilter(display ∘ plotter; seconds=5.0)


stats = TimeFilter(ProgressShower(T), saver, displayer; times)
@time run_RD!(s, M, T; stats, rng)
#savevideo("video2.mp4", saver.stack, plotter)


