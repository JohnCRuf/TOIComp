using Pkg
using LinearAlgebra,DataFrames,DataFramesMeta,Chain, Plots, Cairo, Fontconfig, PlotlyJS, SigFigs
#Same as part b but use delta=0.3
include("problem_2_fns.jl")
#initializing parameters
β = 0.99;
z = 1.5;
α = 0.5;
tol = 1e-6;
δ = 1.0;
k_max = z^(1/(1-α));
k_min = 0.1;
k_grid = collect(range(k_min, k_max, length=500));
v0 = zeros(length(k_grid));
#computing theoretical value function and policy function
v_true = v_theory.(k_grid, β, α, z, δ)
k_true = policy_theory.(z, α, β, k_grid, δ)
k_ss_true = steady_state(z, α, β, δ)[1]
#Using Brute-force algorithm to find the value function

comp_time = @elapsed v_star, k_policy = brute_force(v0, z, β, α, tol, k_grid, δ)
k_ss = computed_steady_state(k_grid, k_policy)
MSE = sum((v_star .- v_true).^2)/length(v_star) 

#Plotting Results and saving outputs
p = Plots.plot(k_grid, v_star, label = ["v_star"], color = ["blue"])
Plots.savefig(p, "../output/part_b_value_comparison_delta_0p3.png")

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
PlotlyJS.savefig(t, "../output/part_b_table.png")
