using Pkg
using LinearAlgebra,DataFrames,DataFramesMeta,Chain, Plots, Cairo, Fontconfig, PlotlyJS, SigFigs
#Same as part b but use delta=0.3
include("problem_3_fns.jl")
#initializing parameters
β = 0.99;
z = 1.5;
α = 0.5;
tol = 1e-6;
δ = 0.3;
k_max = z^(1/(1-α))/δ;
k_min = 0.1;
k_grid = collect(range(k_min, k_max, length=500));
v0L = zeros(length(k_grid));
v0H = zeros(length(k_grid));
#computing theoretical value function and policy function
v_true = v_theory.(k_grid, β, α, z, δ)
k_true = policy_theory.(z, α, β, k_grid, δ)
k_ss_true = steady_state(z, α, β, δ)[1]
#Using Brute-force algorithm to find the value function

comp_time_1 = @elapsed v_starL_1,v_starH_1, k_policyL_1,k_policyH_1 = brute_force(v0L,v0H, z, β, α, tol, k_grid, δ, 1.0)
k_ss_1 = computed_steady_state(k_grid, k_policyH_1)
MSE_1 = sum((v_starH_1 .- v_true).^2)/length(v_starH_1) 

comp_time_0p5 = @elapsed v_starL_0p5,v_starH_0p5, k_policyL_0p5,k_policyH_0p5 = brute_force(v0L,v0H, z, β, α, tol, k_grid, δ, 0.5)
k_ss_0p5 = computed_steady_state(k_grid, k_policyH_0p5)
MSE_0p5 = sum((v_starH_0p5 .- v_true).^2)/length(v_starH_0p5) 

#Plotting Results and saving outputs
p = Plots.plot(k_grid, [v_starH_1, v_starH_0p5], label = ["Gamma = 1 Value Function" "Gamma = 0.5 Value Function"], color = ["blue" "red"],ls = [:dot :solid], legend = :outerbottom, legendcolumns = 2)
Plots.savefig(p, "../output/problem_3_part_c_value_comparison_delta_0p3.png")
