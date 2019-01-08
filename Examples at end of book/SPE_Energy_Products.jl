using JuMP, Clp, Complementarity

producers = 1:4
consumers = 1:5

Supply = [20 30 35 15]
Demand = [20 25 10 10 35]
costs = [10 20 30 5 6;
		 3  15 12 5 6;
		 17 40 5  3 7;
		 50 23 28 12 4]

m = Model(solver = ClpSolver())

@variable(m, F[i in producers, j in consumers] >= 0)
println(1)
@objective(m, Min, sum(costs[i,j]*F[i,j] for i in producers, j in consumers))
@constraint(m, supplyConst[i in producers], Supply[i] - sum(F[i,j] for j in consumers) >= 0)
println(2)
@constraint(m, demandConst[j in consumers], sum(F[i,j] for i in producers) - Demand[j] >= 0)

print(m)
status = solve(m)
println(getvalue(F))
