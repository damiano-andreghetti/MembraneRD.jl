using MembraneRD, Random, Colors, MembraneRD.Filters

species = @species A B EA EB

struct Measure
	time::Float64
    phi_av::Float64
    bc::Float64
    cytoEA::Float64
    cytoEB::Float64
    rho::Float64
end

function Measure(M::Model, s::State, t)
    Nc = nsites(M)
    #calculate phi at each cell (normalized, goes from -1, to 1)
    ϕ(i) = (s.membrane[i,B] - s.membrane[i,A])/(s.membrane[i,A] + s.membrane[i,B] + 1e-10)
    @views totA, totB = sum(s.membrane[:,A]), sum(s.membrane[:,B])
    #measure phi
    ϕav = (totB - totA) / (totA + totB)
    #in some denominators I added + 10^-10 in order to avoid NaNs
    m2 = sum((ϕ(i)-ϕav)^2 for i in 1:Nc) / Nc
    m4 = sum((ϕ(i)-ϕav)^4 for i in 1:Nc) / Nc
    #measure Binder cumulant
    bc = 1 - (m4 / (m2 ^ 2 + 1e-10)) / 3			
    #the following is to check for convergence of rho to 1
    rho = M.rho_0 * ((M.det[1][2] / M.att[1][3]) + totA) / ((M.det[2][2] / M.att[2][3]) + totB)
    Measure(t, ϕav, bc, s.cytosol[EA], s.cytosol[EB], rho)
end

#make measures here
function Measurer(M::Model; name, Nsave)
	measures = Measure[]
    next::Int = 1
    function take_measure(t, s)
        m = Measure(M, s, t)
        push!(measures, m)
        if next % Nsave ==0
            println("T = $t and <ϕ>/c = $(M.ϕav)")
            println("Binder cumulant: $(M.bc)")
            save("$name/config_T=$(round(t,digits=3)).jld",compress=true, "state", s)
            save("$name/measures.jld",compress=true, "measures", measures)
        end
        next += 1
    end
end

function build_model(L)
    G,posx,posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    unit_length = 1.0
    theta = 0.2126944621086619 #Stot/V
    Ac = unit_length^2 #area associated to each cell, set unit lengthscale
    Stot = N*Ac #area of each cell is assumed to be 1, setting unit lengthscale
    V = Stot/theta
    dA = dB = dEA = dEB = 0.01 #this sets unit timescale
    #rates in the theory (do not correspond excatly to those used in the simulations), some slight dimensional changes are needed
    kAa_th = 1.0; kAd_th = 1.0*kAa_th; kAc_th = 1.0
    kBa_th = 1.0; kBd_th = 1.0*kBa_th; kBc_th = 1.3*kAc_th
    KMMA_th = 1.0; KMMB_th = 1.0
    #rates to implement
    kAc = kAc_th; kBc = kBc_th; kAa = kAa_th/V
    kAd = kAd_th; kBa = kBa_th/V; kBd = kBd_th
    KMMA = KMMA_th*Ac; KMMB = KMMB_th*Ac


    cat = ((@catalytic (kAc,KMMA) EA+B => EA+A), 
           (@catalytic (kBc,KMMB) EB+A => EB+B))
    att = ((EA,A,kAa), (EB,B,kBa))
    det = ((EA,kAd),(EB,kBd))
    dif = ((A,dA),(B,dB),(EA,dEA),(EB,dEB))

    M = Model(; species, G, cat, att, det, dif, rho_0 = 0.0)

    totmol = N * 10
    totA, totB = floor(Int, totmol/2), floor(Int, totmol/2)
    totEA, totEB = floor(Int, 0.1*N), floor(Int, 0.1*N)
    memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+theta/(kAd_th/kAa_th)*totA/Stot))
    memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+theta/(kBd_th/kBa_th)*totB/Stot))
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    (; M, mem = [totA, totB, memEA, memEB], cyto = [0.0, 0.0, cytoEA, cytoEB], posx, posy)
end

T = 20000.0
Tmeas = 100.0
Nsave = 10
L = 80

#for reproducibility
seed = 22
rng = Random.Xoshiro(seed)
(; M, mem, cyto, posx, posy) = build_model(L)
s = State(M, mem, cyto; rng)

times = 0:100:T
saver = Pusher(Tuple{Float64,State})
colors = [color("yellow"),color("blue"),color("black"),color("black")]/30
stats = TimeFilter(ProgressShower(T), 
#   Measurer(M; name="test_example_1", Nsave),
#   StopWatchFilter(display ∘ Plotter(posx, posy; colors); seconds=1.0),
    saver; times)
@time run_RD!(s, M, T; stats, rng)

#Measure.(Ref(M), last.(saver.stack), first.(saver.stack))[end]