using Pkg
using LinearAlgebra,DataFrames,DataFramesMeta,Chain, Plots, Cairo, Fontconfig
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
#Using Brute-force algorithm to find the value function
v_star, k_policy = brute_force(v0, z, β, α, tol, k_grid, δ)


#Plotting Results and saving outputs
df = DataFrame(k_grid = k_grid[:], v_star = v_star[:], v_true = v_true[:])
p = Plots.plot(k_grid, [v_star v_true], label = ["v_star" "v_true"], color = ["blue" "red"])
savefig(p, "../output/part_a_value_comparison.png")

df = DataFrame(k_grid = k_grid[:], k_policy = k_policy[:], k_true = k_true[:])
p = Plots.plot(k_grid, [k_policy k_true], label = ["Computed Policy Function" "True Policy Function"], color = ["blue" "red"])
savefig(p, "../output/part_a_policy_comparison.png")
