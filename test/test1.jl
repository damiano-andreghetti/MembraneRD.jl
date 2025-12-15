using MembraneRD, Test, StableRNGs

species = @species A B EA EB

function build_model(L)
    G,layout = MembraneRD.gen_hex_lattice(L)
    N = nv(G)
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
    kAc, kAd, kAa = kAc_th, kAd_th, kAa_th / V
    kBc, kBd, kBa = kBc_th, kBd_th, kBa_th / V
    KMMA, KMMB = KMMA_th*Ac, KMMB_th*Ac
    kAb, kBb = 1e-4, 1e-4


    cat = ((@catalytic (kAc,KMMA) EA+B => EA+A), 
           (@catalytic (kBc,KMMB) EB+A => EB+B))
    att = ((EA,A,kAa), (EB,B,kBa))
    det = ((EA,kAd),(EB,kBd))
    dif = ((A,dA),(B,dB),(EA,dEA),(EB,dEB))
    rea = ((A,),(B,),kAb), ((B,),(A,),kBb)

    M = Model(; species, G, cat, att, det, dif, rea, rho_0 = 0.0)

    # create initial state distribution
    totmol = 10N
    totA, totB = floor(Int, 0.5totmol), floor(Int, 0.5totmol)
    totEA, totEB = floor(Int, 0.1N), floor(Int, 0.1N)
    memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+theta/(kAd_th/kAa_th)*totA/Stot))
    memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+theta/(kBd_th/kBa_th)*totB/Stot))
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    (; M, mem = [totA, totB, memEA, memEB], cyto = [0.0, 0.0, cytoEA, cytoEB], layout)
end



@testset "reproducibility" begin
    T = 2000.0
    L = 100
    rng = StableRNG(22)
    (; M, mem, cyto, layout) = build_model(L)
    s = State(M, mem, cyto; rng)
    run_RD!(s, M, T; rng)
    @test sum(s.membrane[:,1]) == 36934
    @test sum(s.membrane[:,2]) == 63066
end

