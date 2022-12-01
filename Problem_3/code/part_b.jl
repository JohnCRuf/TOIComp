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
γ = 1.0;
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

comp_time = @elapsed v_starL,v_starH, k_policyL,k_policyH = brute_force(v0L,v0H, z, β, α, tol, k_grid, δ, γ)
k_ss = computed_steady_state(k_grid, k_policyH)
MSE = sum((v_starH .- v_true).^2)/length(v_starH) 

#Plotting Results and saving outputs
p = Plots.plot(k_grid, v_starH, label = ["v_star"], color = ["blue"])
Plots.savefig(p, "../output/problem_3_part_b_value_comparison_delta_0p3.png")

#creating table for output
comp_time = round(comp_time, digits=5)
MSE = SigFigs.significantfigures(SigFig(MSE, 5))
k_ss = round(k_ss, digits=5)
k_ss_true = round(k_ss_true, digits=5)
t =PlotlyJS.plot(table(
    header_values = ["Parameter", "Value"],
    cells_values = [
        ["Computation Time", "MSE", "Theoretical Steady State", "Computed Steady State"],
        [comp_time, MSE, k_ss_true, k_ss]
    ]
))
PlotlyJS.savefig(t, "../output/problem_3_part_b_table.png")