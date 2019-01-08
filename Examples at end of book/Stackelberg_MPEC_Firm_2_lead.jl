#Code for the Stackelberg MPEC with Firm 2 as leader
using JuMP, Ipopt, Complementarity


m = Model(solver=IpoptSolver())

@variable(x2())