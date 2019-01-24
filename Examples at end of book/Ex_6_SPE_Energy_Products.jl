#Example 4.4.1.1
using JuMP, Clp

#Here we solve the planners problem as a simple linear program
#the model has n producers and m consumers, and the costs or production are C[n,m]
#Each firm has a capacity constraint given in Supply and each demand region has total that must be met.
#The optimization problem is to meet demand at minimum cost given supply constraints
producers = 1:4
consumers = 1:5

Supply = [20 30 35 15]
Demand = [20 25 10 10 35]
costs = [10 20 30 5 6;
		 3  15 12 5 6;
		 17 40 5  3 7;
		 50 23 28 12 4]

function LP_Model()
	m = Model(solver = ClpSolver())
	@variable(m, F[i in producers, j in consumers] >= 0)
	@objective(m, Min, sum(costs[i,j]*F[i,j] for i in producers, j in consumers))
	@constraint(m, supplyConst[i in producers], Supply[i] - sum(F[i,j] for j in consumers) >= 0)
	@constraint(m, demandConst[j in consumers], sum(F[i,j] for i in producers) - Demand[j] >= 0)

	# print(m)
	status = solve(m)
	println(getvalue(F))
end

LP_Model()
