# PlatPos
Simulation for bidder's strategy which contains solving HJB equation, transportation equation and fix poiont problem.

Our simulation can be divided into four parts as following,
First step(Z_eq1.m) Sovle Z_eq by Simpson's Rule which is trivial;

Second step(HJB1.m) Sovle HJB equation by finite difference method and get strategy v;

Third (v_star.m) Trivial;

Fourth (M_upwind1.m) We use upwind scheme to solve the transport equation to get density M(t,x);

Fifth (Int.m)Finally we sovle the fix point equation and finish the first loop.

Then perform above step to get convergence.

The main code is run.m


#####################################################
#interpretation
...
