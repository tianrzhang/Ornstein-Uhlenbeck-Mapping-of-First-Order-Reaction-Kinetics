import numpy as np
import matplotlib.pyplot as plt

# Mapping O-U simulation onto a first-order reaction kinetics (α-D-glucose ⇌ β-D-glucose)
# O-U equation is as follows: dX_t = θ(μ-X_t) + σdW_t
# All literature values are pulled from https://pubs.acs.org/doi/abs/10.1021/jp906523s
# Note that literature value of k_1 = k_backward (β to α) and k_2 = k_forward (α to β) in this project
# The experiment bearing the literature values was conducted in tridistilled water
# Instead of mapping concentration, we can map type of molecule count instead i.e. number of α-D-glucose molecules

k_forward = 4.91*10**-5
k_backward = 2.76*10**-5
k_total = k_forward + k_backward
K_i = k_forward / k_backward
p = 1/(1 + K_i)
theta = k_total

def mean_revert(N, alpha_count, dt):
    alpha_initial = N
    miu = (k_backward * alpha_initial)/theta
    return theta*(miu - alpha_count)*dt

def noise(N, dt):
    sigma = np.sqrt(2*theta*N*p*(1-p))
    return sigma*np.sqrt(dt)*np.random.normal(0,1)

def reaction_fixed(N, alpha_count, t_max, dt):
    """Here, alpha_count is the number of α-D-glucose molecules
    at time = 0"""
    steps_count = int((t_max/dt)+1)
    data = np.zeros(steps_count)
    data[0] = alpha_count
    for i in range(1, steps_count):
        alpha_count = alpha_count + mean_revert(N, alpha_count, dt) + noise(N, dt)
        data[i] = alpha_count
    return data

def reaction_hit(N, alpha_count, dt, limit, expire):
    """Here, alpha_count is the number of α-D-glucose molecules
        at time = 0 as well"""
    reached_limit = ""
    time = 0
    data = []
    constraint = expire * 24 * 60 * 60
    while alpha_count > limit and time < constraint:
        time += dt
        alpha_count = alpha_count + mean_revert(N, alpha_count, dt) + noise(N, dt)
        data.append(alpha_count)
    if alpha_count <= limit:
        reached_limit = "reaction has reached desired α-D-glucose molecule count"
    else:
        reached_limit = "the desired α-D-glucose molecule count was not reached within the designated time frame"
    return data, reached_limit, time

def simulate_ou_fixed(N, t_max, dt):
    alpha_count = N
    data = reaction_fixed(N, alpha_count, t_max, dt)
    return data


def simulate_ou_hit(N, dt, limit, expire):
    alpha_count = N
    data, reached_limit, time = reaction_hit(N, alpha_count, dt, limit, expire)
    return data, reached_limit, time

def plotting_sim_fixed(N, t_max, dt):
    plt.plot(range(int(t_max/dt)+1), simulate_ou_fixed(N, t_max, dt), color = "cornflowerblue")
    plt.xlabel("Steps taken (time/dt)")
    plt.ylabel("α-D-glucose molecule count")
    plt.title("α-D-glucose molecule count against time")
    plt.show()
    return

plotting_sim_fixed(10, 1000000, 50)
plotting_sim_fixed(100, 1000000, 50)
plotting_sim_fixed(1000, 1000000, 50)
plotting_sim_fixed(10000, 1000000, 50)

def plotting_sim_hit(N, limit, dt, expire):
    data, reached_limit, time = simulate_ou_hit(N, dt, limit, expire)
    plt.plot(range(len(data)), data, color = "cornflowerblue")
    plt.xlabel("Steps taken (time/dt)")
    plt.ylabel("α-D-glucose molecule count")
    plt.title("Time taken for α-D-glucose molecule count to reach desired count")
    plt.show()
    print(reached_limit)
    return

plotting_sim_hit(10, 2, 50, 7)
plotting_sim_hit(100, 20, 50, 7)
plotting_sim_hit(1000, 200, 50, 7)
plotting_sim_hit(10000, 2000, 50, 7)

# We seek to validate this phenomenon we observe from the simulation above: relative noise decreases as N increases and absolute noise increases with N
# This can be done quantitively via finding the difference between the noisy simulation and noiseless exponential decay i.e. intended shape of the process
# Why are we using the noiseless exponential decay? This is because it maps out the desired path of the process: relaxing towards the equilibrium
# Therefore any fluctuations around this intended path can be classified as "Noise"
# Given our reaction/process is binomial, we can declare that the spread around out intended path values (standard deviation) = sqrt(Np(1-p))
# Intended outcome is to yield two graphs, one of absolute and one of relative standard deviation and find the noise by finding the difference between the values (residuals)

def exp_path(N, t_max, dt):
    miu = N*p
    steps_count = int((t_max/dt)+1)
    time_term = np.arange(0, steps_count, 1) * dt
    alpha_decaying = N*np.exp(-theta*time_term) + miu*(1 - np.exp(-theta*time_term))
    return alpha_decaying

def plotting_exp_path(N, t_max, dt):
    alpha_decaying = exp_path(N, t_max, dt)
    steps_count = int((t_max/dt+1))
    steps = np.arange(0, steps_count, 1)
    plt.plot(steps, alpha_decaying, color = "cornflowerblue")
    plt.xlabel("Steps taken (time/dt)")
    plt.ylabel("α-D-glucose molecule count")
    plt.title("Exponential decay of α-D-glucose molecule count to equilibrium count")
    plt.show()
    return

plotting_exp_path(10, 100000, 50)
plotting_exp_path(100, 100000, 50)
plotting_exp_path(1000, 100000, 50)
plotting_exp_path(10000, 100000, 50)

def plotting_compare(N, t_max, dt):
    steps_count=int((t_max/dt)+1)
    alpha_decaying = exp_path(N, t_max, dt)
    steps = np.arange(0, steps_count, 1)
    plt.plot(steps, alpha_decaying, label = "exponential path", color = "cornflowerblue", alpha = 0.5)
    fixed_sim = simulate_ou_fixed(N, t_max, dt)
    plt.plot(steps, fixed_sim, label = "O-U simulation path", color = "firebrick", alpha = 0.5)
    plt.xlabel("Steps taken (time/dt)")
    plt.ylabel("α-D-glucose molecule count")
    plt.title("Exponential decay and O-U simulation of α-D-glucose molecule count reaching equilibrium", fontsize = 9)
    plt.legend()
    plt.show()
    return

plotting_compare(10, 100000, 50)
plotting_compare(100, 100000, 50)
plotting_compare(1000, 100000, 50)
plotting_compare(10000, 100000, 50)

def residuals(N, t_max, dt):
    alpha_decaying = exp_path(N, t_max, dt)
    fixed_sim = simulate_ou_fixed(N, t_max, dt)
    fluctuation = fixed_sim - alpha_decaying
    return fluctuation

def plotting_residuals(N, t_max, dt):
    steps_count=int((t_max/dt)+1)
    steps = np.arange(0, steps_count, 1)
    fluctuation = residuals(N, t_max, dt)
    plt.plot(steps, fluctuation, color = "cornflowerblue")
    plt.xlabel("Steps taken (time/dt)")
    plt.ylabel("Fluctuations around the deterministic path")
    plt.title("Fluctuations around the deterministic path against steps taken")
    plt.show()
    return

plotting_residuals(10, 100000, 50)
plotting_residuals(100, 100000, 50)
plotting_residuals(1000, 100000, 50)
plotting_residuals(10000, 100000, 50)

# In order for us to know whether the model (exponential path) we are using is unbiased or not, we must compare the mean of the fluctuations against SDOM
# SDOM is given by σ/sqrt(N) } this is however for data with independent points
# Since our data points are all affecting one another, i.e. they are correlated, we must use SODM = σ/sqrt(ESS), whereby ESS is the number of independent points carrying the same information as one correlated sample
# ESS is given by duration of process divided by 2 * correlation time (τ_c)
# correlation time (τ_c) is given by comparing the general definition of ACF exp(-τ/τ_c) and the OU specific ACF exp(-k*abs(τ))
# This gives us -k*abs(τ) = -τ/τ_c, resulting in k = 1/τ_c whereby k is theta in our case

def scatter_percent(percentage, target):
    """ This function computes out the duration time needed and the effective sample size for different targets,
    standard deviation (std) or mean via the percentage spread desired """
    tau_c = 1/theta
    decimal = percentage/100
    if target == "mean":
        eff_samp = 1/(decimal)**2
    elif target in ("std", "standard deviation"):
        eff_samp = 1/(2*(decimal)**2)
    t_max = eff_samp*(2*tau_c)
    return eff_samp, t_max

# For this project, we want a 1% scatter percentage for both std and mean

print(scatter_percent(1, "std"))

def noise_against_N(N_vals, t_max, dt):
    abs_std_vals = []
    rel_std_vals = []
    for N in N_vals:
        fluctuations = residuals(N, t_max, dt)
        abs_std_vals = abs_std_vals + [np.std(fluctuations)]
        means = N*p
        rel_std_vals = rel_std_vals + [np.std(fluctuations)/means]
    return abs_std_vals, rel_std_vals

print(noise_against_N([10,100,1000,10000], scatter_percent(1, "std")[1], 500))

def plotting_log_abs(N_vals, dt):
    t_max = scatter_percent(1, "std")[1]
    std_vals = np.log10(noise_against_N(N_vals, t_max, dt)[0])
    plt.plot(np.log10(N_vals), std_vals, color = "cornflowerblue")
    plt.xlabel("Log(N)")
    plt.ylabel("Log(standard deviation of residuals)")
    plt.title("Log(standard deviation) against Log(N)")
    plt.show()
    return

plotting_log_abs([10,100,1000,10000], 500)

def plotting_log_rel(N_vals, dt):
    t_max = scatter_percent(1, "std")[1]
    std_vals = np.log10(noise_against_N(N_vals, t_max, dt)[1])
    plt.plot(np.log10(N_vals), std_vals, color = "cornflowerblue")
    plt.xlabel("Log(N)")
    plt.ylabel("Log(standard deviation of residuals/mean)")
    plt.title("Log(standard deviation/mean) against Log(N)")
    plt.show()
    return

plotting_log_rel([10,100,1000,10000], 500)

print(np.polyfit(np.log10([10,100,1000,10000]), np.log10(noise_against_N([10,100,1000,10000], scatter_percent(1, "std")[1], 500)[0]), 1))

print(np.polyfit(np.log10([10,100,1000,10000]), np.log10(noise_against_N([10,100,1000,10000], scatter_percent(1, "std")[1], 500)[1]), 1))

def standard_deviation(N_vals):
    absolute = []
    relative = []
    for N in N_vals:
        std_abs = np.sqrt(N*p*(1-p))
        std_rel = np.sqrt(N*p*(1-p))/(N*p)
        absolute = absolute + [std_abs]
        relative = relative + [std_rel]
    return absolute, relative

def compare_abs_plot(N_vals, dt):
    t_max = scatter_percent(1, "std")[1]
    std_vals = np.log10(noise_against_N(N_vals, t_max, dt)[0])
    std_abs = np.log10(standard_deviation(N_vals)[0])
    plt.plot(np.log10(N_vals), std_vals, label = "measured", color = "cornflowerblue", alpha = 0.5)
    plt.plot(np.log10(N_vals), std_abs, label = "derived", color = "firebrick", alpha = 0.5)
    plt.xlabel("Log(N)")
    plt.ylabel("Log(standard deviation of residuals)")
    plt.title("Comparison between derived and measured Log(std) against Log(N)")
    plt.legend()
    plt.show()
    return

compare_abs_plot([10,100,1000,10000], 500)

def compare_rel_plot(N_vals, dt):
    t_max = scatter_percent(1, "std")[1]
    std_vals = np.log10(noise_against_N(N_vals, t_max, dt)[1])
    std_rel = np.log10(standard_deviation(N_vals)[1])
    plt.plot(np.log10(N_vals), std_vals, label = "measured", color = "cornflowerblue", alpha = 0.5)
    plt.plot(np.log10(N_vals), std_rel, label = "derived", color = "firebrick", alpha = 0.5)
    plt.xlabel("Log(N)")
    plt.ylabel("Log(standard deviation of residuals/mean)")
    plt.title("Comparison between derived and measured Log(std/mean) against Log(N)", fontsize = 10)
    plt.legend()
    plt.show()
    return

compare_rel_plot([10,100,1000,10000],500)

# In quant, autoregressive (AR) regression is used to predict future behaviour of a model from historical data
# Since this project still has the background aim to link chemistry to finance, we will be recovering parameters via AR regression and predicting future behaviour via the recovery
# The function we wish to build now recovers θ, μ (both from our mean reversion term), and σ (from the nosie term) for our OU model
# Embedded in the slope and the intercept are the values of θ and μ
# The residual data given by the difference of the actual path and predicted path can be used to derive σ
# However, given we only have data generated from known parameters, we can only use those to "simulate" parameter recovery

def recover_params(data, dt):
    model_predict = np.polyfit(data[:-1], data[1:], 1)
    m = model_predict[0]
    c = model_predict[1]
    θ = (1-m)/dt
    μ = c/(1-m)
    path_predict = m*data[:-1] + c
    actual_path = data[1:]
    residuals = actual_path - path_predict
    σ = np.std(residuals)/np.sqrt(dt)
    return θ, μ, σ

data_10_N = reaction_fixed(10, 10, scatter_percent(1, "std")[1], 500)
data_100_N = reaction_fixed(100, 100, scatter_percent(1, "std")[1], 500)
data_1000_N = reaction_fixed(1000, 1000, scatter_percent(1, "std")[1], 500)
data_10000_N = reaction_fixed(10000, 10000, scatter_percent(1, "std")[1], 500)

print("The following data gives θ, μ, and σ in that order")

print(recover_params(data_10_N, 500))
print(recover_params(data_100_N, 500))
print(recover_params(data_1000_N, 500))
print(recover_params(data_10000_N, 500))

# We can further validate the hypothesis that at large N counts, noise dies down by comparing the results given by the OU process and Gillespie
# The Gillespie function from a previous project is reused here, along with the functions cutting off the transient portions

def Gillespie(N, t_max):
    t = 0
    alpha_count = N
    beta_count = 0
    time = []
    alpha_mols = []
    while t < t_max:
        a1 = k_forward * alpha_count
        a2 = k_backward * beta_count
        a0 = a1 + a2
        if a0 == 0:
            break
        τ = np.random.exponential(1/a0)
        choose = np.random.uniform(0, 1)
        if choose < a1/a0:
            alpha_count -= 1
            beta_count += 1
        else:
            alpha_count += 1
            beta_count -= 1
        t += τ
        time.append(t)
        alpha_mols.append(alpha_count)
    return time, alpha_mols

def stationary_count_split(data, N, window, flat, repeat):
    previous_mean = 0
    current_mean = 0
    persist_count = 0
    for i in range(2*window, len(data)):
        current_mean = sum(data[i-window:i]) / window
        previous_mean = sum(data[i-2*window:i-window]) / window
        derivative = current_mean - previous_mean
        if abs(derivative/N) < flat:
            persist_count += 1
        else:
            persist_count = 0
        if persist_count >= repeat:
            return i
    return len(data) // 2

def stationary_count(N, t_max, window, flat, repeat):
    time, alpha_mols = Gillespie(N, t_max)
    split = stationary_count_split(alpha_mols, N, window, flat, repeat)
    return time[split:], alpha_mols[split:]

stationary_10_gillespie = stationary_count(10, 130000, 50, 0.5, 5)
stationary_100_gillespie = stationary_count(100, 130000, 50, 0.01, 10)
stationary_1000_gillespie = stationary_count(1000, 130000, 50, 0.001, 25)
stationary_10000_gillespie = stationary_count(10000, 130000, 50, 0.00001, 50)

mols_10_fractions = [mol/10 for mol in stationary_10_gillespie[1]]
mols_100_fractions = [mol/100 for mol in stationary_100_gillespie[1]]
mols_1000_fractions = [mol/1000 for mol in stationary_1000_gillespie[1]]
mols_10000_fractions = [mol/10000 for mol in stationary_10000_gillespie[1]]

shared_edges = np.arange(0.2, 0.6, 0.005)

plt.hist(mols_10_fractions, bins = shared_edges, density = True, label = "N=10", alpha = 0.5)
plt.hist(mols_100_fractions, bins = shared_edges, density = True, label = "N=100", alpha = 0.5)
plt.hist(mols_1000_fractions, bins = shared_edges, density = True, label = "N=1000", alpha = 0.5)
plt.hist(mols_10000_fractions, bins = shared_edges, density = True, label = "N=10000", alpha = 0.5)
plt.xlabel("α-D-glucose molecule in proportion to N")
plt.ylabel("Probability Density")
plt.title("A distribution of α-D-glucose molecules in proportion to N at different N (Gillespie)", fontsize = 10)
plt.legend()
plt.show()

# For the OU process, using the stationary_count_split function is less appropriate given the fluctuations around the intended path
# This means we should just manually cut off the data for the OU process, leaving only the data @ equilibrium (stationary) point
# Relaxation time, or rather time for the process to forget its starting condition is given by (1/θ)/dt

def ou_data_split(data, dt):
    cutoff = int(5*(1/theta)/dt)
    split_data = data[cutoff:]
    return split_data

split_10 = ou_data_split(data_10_N, 500)
split_100 = ou_data_split(data_100_N, 500)
split_1000 = ou_data_split(data_1000_N, 500)
split_10000 = ou_data_split(data_10000_N, 500)

split_10_fractional = [mol/10 for mol in split_10]
split_100_fractional = [mol/100 for mol in split_100]
split_1000_fractional = [mol/1000 for mol in split_1000]
split_10000_fractional = [mol/10000 for mol in split_10000]

# A comparison between Gillespie and OU distributions of α-D-glucose molecules at equilibrium

shared_edges_2 = np.arange(0, 0.8, 0.005)
fig, axs = plt.subplots(2,2)

axs[0,0].hist(split_10_fractional, bins = shared_edges_2, density = True, label = "Ornstein-Uhlenbeck", alpha = 0.5, color = "firebrick")
axs[0,0].hist(mols_10_fractions, bins = shared_edges_2, density = True, label = "Gillespie", alpha = 0.5, color = "cornflowerblue")
axs[0,0].set_title("N=10", fontsize = 9)
axs[0,0].legend(fontsize = "xx-small", loc = "upper right")

axs[0,1].hist(split_100_fractional, bins = shared_edges_2, density = True, label = "Ornstein-Uhlenbeck", alpha = 0.5, color = "firebrick")
axs[0,1].hist(mols_100_fractions, bins = shared_edges_2, density = True, label = "Gillespie", alpha = 0.5, color = "cornflowerblue")
axs[0,1].set_title("N=100", fontsize = 9)
axs[0,1].legend(fontsize = "xx-small", loc = "upper right")

axs[1,0].hist(split_1000_fractional, bins = shared_edges_2, density = True, label = "Ornstein-Uhlenbeck", alpha = 0.5, color = "firebrick")
axs[1,0].hist(mols_1000_fractions, bins = shared_edges_2, density = True, label = "Gillespie", alpha = 0.5, color = "cornflowerblue")
axs[1,0].set_title("N=1000", fontsize = 9)
axs[1,0].legend(fontsize = "xx-small", loc = "upper right")

axs[1,1].hist(split_10000_fractional, bins = shared_edges_2, density = True, label = "Ornstein-Uhlenbeck", alpha = 0.5, color = "firebrick")
axs[1,1].hist(mols_10000_fractions, bins = shared_edges_2, density = True, label = "Gillespie", alpha = 0.5, color = "cornflowerblue")
axs[1,1].set_title("N=10000", fontsize = 9)
axs[1,1].legend(fontsize = "xx-small", loc = "upper right")

fig.suptitle("Comparison between distribution of α-D-glucose molecules in Gillespie vs OU process at different N", fontsize = 8.25)
fig.supxlabel("α-D-glucose molecule in proportion to N", fontsize = 9)
fig.supylabel("Probability Density", fontsize = 9)
plt.tight_layout()
plt.show()