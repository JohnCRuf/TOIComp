using Pkg
using LinearAlgebra,DataFrames,DataFramesMeta,Chain, Plots, Cairo, Fontconfig
#same as part b but now change grid size and track computational time
using Pkg
using LinearAlgebra,DataFrames,DataFramesMeta,Chain, Plots, Cairo, Fontconfig, PlotlyJS, SigFigs
#Same as part b but use delta=0.3
include("problem_2_fns.jl")
#initializing parameters
β = 0.99;
z = 1.5;
α = 0.5;
tol = 1e-6;
δ = 0.3;
k_max = z^(1/(1-α))/δ;
k_min = 0.1;
g = [50, 100, 500, 1000, 2000] #grid sizes
#Using Brute-force algorithm to find the value function
comp_time, err = comp_and_error(k_min, k_max, g, z, α, β, δ, tol)

#Plotting Results and saving outputs
p = Plots.plot(g, comp_time, label = ["Computing Time"], color = ["black"])
Plots.savefig(p, "../output/part_c_computation_time.png")

p2 = Plots.plot(g, err, label = ["Steady State Error"], color = ["green"])
Plots.savefig(p2, "../output/part_c_error.png")

