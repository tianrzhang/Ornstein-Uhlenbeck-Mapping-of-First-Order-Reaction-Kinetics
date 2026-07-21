# Ornstein-Uhlenbeck-Mapping-of-First-Order-Reaction-Kinetics

Stochastic model of glucose mutarotation, validated against deterministic rate law and stochastic kinetics of Gillespie process. 

This project utilises the mean-reverting, Gauss-Markov nature of the Ornstein-Uhlenbeck process to simulate the reaction kinetics of the first order equilibrium reaction of Glucose mutarotation. It was demonstrated that the Ornstein-Uhlenbeck process reproduces both the deterministic rate law trajectory and the stochastic kinetics of the Gillespie approach at macroscopic levels.

To numerically process the Ornstein-Uhlenbeck model, the Euler-Maruyama scheme was used so that the SDE of the process could be discretised. Key parameters that could be used globally were mapped out first, followed by individual terms of the function, mean-reversion and noise. 

