using FractionalKrylovMC, Plots, Statistics

trajectories = 5
nsimvec = Int.([2^i for i in 10:21]);
experiments = length(nsimvec)
k = 25 # size of matrices
cutoff = 1e8
err_to_plot = zeros(experiments,trajectories)

for trajectory in 1:trajectories
    @show trajectory
    α = rand()*0.1 + 0.9; γ = rand()*0.1 + 0.9
    problem = create_random_problem(α,γ,k,MittagLefflerRand())
    for exper in 1:experiments
        MCKrySol = solve(problem,nsimvec[exper],cutoff,MCSolver())
        err_to_plot[exper,trajectory] = relative_error(problem.uT_spectral,MCKrySol.uT)
    end
end

err_to_plot
plot(log.(nsimvec),log.(mean(err_to_plot,dims=2)),title="$(trajectories) trajectories",xlabel="Log sim",ylabel="Mean log relative err")
savefig("test/figs2/moresim_lesserr_plot2.png")



# 1st
# trajectories = 6
# nsimvec = Int.([1e3;1e4;1e5;1e6;1e7]); experiments = length(nsimvec)
# err_to_plot = zeros(experiments,trajectories)
# k = 40 # size of matrices
# cutoff = 1e5
# α = rand()*0.9 + 0.1; γ = rand()*0.9 + 0.1

# 2nd
# trajectories = 20
# nsimvec = Int.([1e1;1e2;1e3;1e4;1e5]); experiments = length(nsimvec)
# err_to_plot = zeros(experiments,trajectories)
# k = 30 # size of matrices
# cutoff = 1e7
# α = rand()*0.4 + 0.6; γ = rand()*0.4 + 0.6

# 3rd
# trajectories = 30
# nsimvec = Int.([2^i for i in 5:14]); experiments = length(nsimvec)
# err_to_plot = zeros(experiments,trajectories)
# k = 25 # size of matrices
# cutoff = 1e7

# 4th
# trajectories = 20
# nsimvec = Int.([2^i for i in 5:20]); experiments = length(nsimvec)
# k = 15 # size of matrices
# cutoff = 1e7

# 5th
# trajectories = 20
# nsimvec = Int.([2^i for i in 14:21]); experiments = length(nsimvec)
# k = 25 # size of matrices
# cutoff = 1e8
# α = rand()*0.1 + 0.9; γ = rand()*0.1 + 0.9
