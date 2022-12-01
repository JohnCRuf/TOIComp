function bellman(v_L, v_H, k_grid, z, β, α, δ, γ)
    if length(v_H) != length(k_grid) | length(v_L) != length(k_grid)
        error("v and k_grid must be the same length")
    end
    v_new_L = zeros(length(k_grid))
    v_new_H = zeros(length(k_grid))
    k_chosen_L = zeros(length(k_grid))
    k_chosen_H = zeros(length(k_grid))
    for i in 1:length(k_grid)
            k = k_grid[i]
            obj_L = [k_grid[j] <= 0.5*k^α+ (1-δ)*k ? log(0.5*k^α + (1-δ)*k- k_grid[j]) + β*(γ*v_H[j] + (1-γ)*v_L[j]) : -Inf for j = 1:length(k_grid)]
            obj_H = [k_grid[j] <= 1.5*k^α+ (1-δ)*k ? log(1.5*k^α + (1-δ)*k- k_grid[j]) + β*(γ* v_H[j] + (1-γ)*v_L[j]) : -Inf for j = 1:length(k_grid)]
            v_new_H[i] = maximum(obj_H)
            v_new_L[i] = maximum(obj_L)
            k_chosen_L[i] = k_grid[findmax(obj_L)[2]] 
            k_chosen_H[i] = k_grid[findmax(obj_H)[2]] 
    end
    return (v_new_L, v_new_H, k_chosen_L, k_chosen_H)
end
#Iterating to Converge to True Value Function
function brute_force(v0L,v0H, z, β, α, tol, k_grid, δ,γ)
err = Inf
while err > tol
    v1L,v1H, k0L, k0H = bellman(v0L,v0H, k_grid, z, β, α, δ, γ)
    err_L = maximum(abs.(v1L .- v0L))
    err_H = maximum(abs.(v1H .- v0H))
    err = max(err_L, err_H)
    v0L = copy(v1L)
    v0H = copy(v1H)
end
v1L,v1H, k0L, k0H = bellman(v0L,v0H, k_grid, z, β, α, δ, γ)
return (v1L,v1H, k0L, k0H )
end
#Defining Theoretical Value Function
function v_theory(k, β, α, z, δ)
    kss, B, A=steady_state(z, α, β, δ)
    return A + B*log(k)
end

#Creating Theoretical Policy Function
function policy_theory(z, α, β, k, δ)
    kss, B, A =steady_state(z, α, β, δ)
    return z*(1-B^(-1)*α)*k^α +(1-B^(-1))*(1-δ)*k
end
#Creating steady State fn
function steady_state(z, α, β, δ)
    kss = ((β^(-1)-(1-δ))*z^(-1)*α^(-1))^(1/(α-1))
    B = (-z*α*kss^(α-1)-(1-δ))/(1-z*kss^(α-1)-(1-δ))
    A = log(z*kss^(α)+ (1-δ)*kss-kss)/(1-β) - B*log(kss) 
return kss, B, A
end 

function computed_steady_state(k_grid, k_policy) 
err_old = Inf
k_ss = 0
    for i in 1:length(k_grid)
        err = abs(k_policy[i]-k_grid[i])
        if err < err_old
            k_ss = k_grid[i]
            err_old = err
        end
    end
return k_ss
end

function comp_and_error(k_min, k_max, grid_sizes, z, α, β, δ, tol)
    i = 1
    comp_time = zeros(length(grid_sizes))
    err = zeros(length(grid_sizes))
    for g in grid_sizes
        k_grid = collect(range(k_min, k_max, length=g));
        v0 = zeros(g);
        c_iter = @elapsed v_star, k_policy = brute_force(v0, z, β, α, tol, k_grid, δ)
        comp_time[i] = c_iter
        k_ss = computed_steady_state(k_grid, k_policy)
        k_ss_true = steady_state(z, α, β, δ)[1]
        err[i] = abs(k_ss - k_ss_true)
        println(err[i])
        i += 1
    end
    return comp_time, err
end