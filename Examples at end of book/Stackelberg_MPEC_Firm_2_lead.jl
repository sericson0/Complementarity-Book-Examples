#Code for the Stackelberg MPEC with Firm 2 as leader
using JuMP, Ipopt, Complementarity

followers = 1:2
production_process = 1:2
C_l = [0.62 1.25]
C_f = [0.55 0.81;
	   0.78 1.35]
K_l = [17000 22000]
K_f = [21000 16000;
       18000 14000]

alpha = 2.5
beta = 0.0000166666667

m = Model(solver=IpoptSolver())

@variable(m, x_l[b in production_process] >= 0)
@variable(m, x_f[i in followers, b in production_process] >= 0)
@variable(m, q >= 0)
@variable(m, gamma[i in followers, b in production_process] >= 0)
##
@objective(m, Min, sum(C_l[b]*x_l[b] - (alpha - beta*q)*x_l[b] for b in production_process))
@constraint(m, q - sum(x_f[i,b] for i in followers, b in production_process) - sum(x_l[b] for b in production_process) == 0) #Markets clear
for b in production_process
	@constraint(m, K_l[b] - x_l[b] >= 0)
end
##
for i in followers, b in production_process
	@complements(m, 0 <= K_f[i,b] - x_f[i,b], gamma[i, b] >= 0)
	@complements(m, 0 <= C_f[i, b] - (alpha - beta*q) + gamma[i, b], x_f[i, b] >= 0)
end

# print(m)
solve(m)

@show getobjectivevalue(m)
Q = getvalue(q)
X = getvalue(x_l)
Price = alpha - beta*Q
ConPmt = Price*Q
ConSurp = alpha*Q - 0.5*beta*Q^2 - ConPmt

Rev_l = Price*sum(X)
VarCost_l = sum(C_l[b]*X[b] for b in production_process)
ProdSurp = Rev_l - VarCost_l

println("Price: ", Price, " leader production: ", X, ", consumer surplus: ", ConSurp, ", leader revenue: ", Rev_l, ", leader surplus: ", ProdSurp)