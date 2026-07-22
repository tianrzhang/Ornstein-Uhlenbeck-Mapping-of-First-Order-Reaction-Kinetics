# Ornstein-Uhlenbeck Mapping of First Order Reaction Kinetics
Stochastic model of glucose mutarotation, validated against deterministic rate law and stochastic kinetics of Gillespie process. 
# Overview
This project utilises the mean-reverting, Gauss-Markov nature of the Ornstein-Uhlenbeck process to simulate the reaction kinetics of the first order equilibrium reaction of Glucose mutarotation. It was demonstrated that the Ornstein-Uhlenbeck process reproduces both the deterministic rate law trajectory and the stochastic kinetics of the Gillespie approach at macroscopic levels.
# Process
To numerically process the Ornstein-Uhlenbeck model, the Euler-Maruyama scheme was used so that the SDE of the process could be discretised. Key parameters that could be used globally were mapped out first, followed by individual terms of the function, mean-reversion and noise. AR(1) regression was then utilised to conduct parameter recovery of core values: theta, miu, sigma.
# Results
![A graph demonstrating at large N values, the deterministic approximation matches the Ornstein-Uhlenbeck model](figure/deterministic_overlay.png)

![A graph demonstrating the relationship between the absolute standard deviation and N on log-log axes](figure/absolute.png)

![A graph demonstrating the relationship between the relative staandard deviation and N on log-log axes](figure/relative.png)

![A figure showing the convergence of the Ornstein-Uhlenbeck to the exact Gillespie model](figure/Gillespie_OU.png)

# Key Findings
Absolute standard deviation increases along with N. For a log-log relationship, for every N, standard deviation increases by 1/2. Absolute standard deviation is proportional to the square root N.
Relative standard deviation decreases along with N. For a log-log relationship, for every N, standard deviation decreases by 1/2. Relative standard deviation is inversely proportional to the square root of N.
At large N values, the macroscopic level, the Ornstein-Uhlenbeck model demonstrates convergence to the Gillespie and matching the deterministic rate law.
# Running the Code
pip install -r requirements.txt
python mutarotation_ou.py
# References
Rate constants: DOI 10.1021/jp906523s

Gardiner, Handbook of Stochastic Methods (OU process, §3.7-3.8)
