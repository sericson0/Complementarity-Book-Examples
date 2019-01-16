using JuMP, Clp, LinearAlgebra
#This is a very simple capacity expansion model with a single node and two production types
Y = 1:10
P = 1:5
hours_in_period = [8760 / length(P) for p in P]
##
D0 = 10000
demand_growth_rate = 0.02
Dbase = [D0*(1+demand_growth_rate)^(y-1) for y in Y]
D_p_mult = [1 + (p - 1)/length(P) for p in P]
Demand = [Dbase[y]*D_p_mult[p] for y in Y, p in P]
##



technologies = ["Baseload" "Peakers"]
S = 1:length(technologies)
var_cost = [20 60]
cap_cost = [1500 1000].*1000
capacity_factors = [1 for s in S]
discount_rate = 0.03
beta = 1 - discount_rate
# ##

m = Model(solver = ClpSolver())

@variable(m, x[s in S, y in Y] >= 0) #Investment in each year for each generation source 
@variable(m, g[s in S, y in Y, p in P] >= 0) #Generation in each year for each generation source for each period
@expression(m, K[s in S, y in Y], sum(x[s, y1] for y1 in 1:y))

@objective(m, Min, sum(beta^(y-1)*(sum(x[s, y] + sum(g[s, y, p] for p in P) for s in S)) for y in Y))
@constraint(m, MeetDemand[y in Y, p in P], sum(g[s, y, p] for s in S) >= Demand[y, p]*hours_in_period[p])
@constraint(m, CapacityConstraint[y in Y, s in S, p in P], K[s, y]*hours_in_period[p]*capacity_factors[s] >= g[s, y, p])
@constraint(m, KConstraint[s in S, y in Y, p in P], K[s, y]*hours_in_period[p]*capacity_factors[s] >= g[s,y,p])

println("start to solve")
solve(m)
##
X = getvalue(x)
G = getvalue(g)
K = [sum(X[s, y1] for y1 in 1:y) for s in S, y in Y]
##
println("Total Cost: ", getobjectivevalue(m))
for y in Y
	println("Year: ", y, " investment: ", X[:,y], " generation: ", sum(G[:, y, p] for p in P), " capacity_factors: ", [sum(G[s,y,p] for p in P)./(K[s,y].*8760) for s in S])
end


