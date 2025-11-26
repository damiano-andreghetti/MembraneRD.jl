using MembraneRD, Random, MembraneRD.Filters, Colors

species = @species A B C EAB EAC EBA EBC ECA ECB


function build_model3(L)
    G, posx, posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    cat = (EBA,B,A,1.0,1.0), (EAB,A,B,1.0,1.0), (EAC,A,C,1.0,1.0), 
          (ECA,C,A,1.0,1.0), (EBC,B,C,1.0,1.0),(ECB,C,B,1.0,1.0)
    att = (EBA,A,1.0), (EBC,C,1.0), (EAB,B,1.0), (EAC,C,1.0), (ECA,A,1.0), (ECB,B,1.0)
    det = (EBA,1.0), (EBC,1.0), (EAB,1.0), (EAC,1.0), (ECA,1.0), (ECB,1.0)
    dif = (((i,0.1) for i in species)...,)
    M = Model(; species, G, cat, att, det, dif, rho_0 = 0.0)
    (; M, mem = [fill(10*N, 3); fill(0, 6)], 
        cyto = [fill(0,3); fill(floor(Int, 0.5*N),6)], 
        posx, posy)
end

T = 5000.0
Nsave = 10
L = 100

rng = Random.Xoshiro(22)
(; M, mem, cyto, posx, posy) = build_model3(L)
s = State(M, mem, cyto; rng)

colors=[RGB(m==1,m==2,m==3)/30 for m in species]
prog = ProgressShower(T)
plot = StopWatchFilter(display ∘ Plotter(posx, posy; colors); seconds=1.0)
stats = TimeFilter(prog, plot; times=0:T/100:T)

@profview run_RD!(s, M, T; stats, rng)


