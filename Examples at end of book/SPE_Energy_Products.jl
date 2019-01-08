#Example 4.4.1.1
using JuMP, Clp, Complementarity

producers = 1:4
consumers = 1:5

Supply = [20 30 35 15]
Demand = [20 25 10 10 35]
costs = [10 20 30 5 6;
		 3  15 12 5 6;
		 17 40 5  3 7;
		 50 23 28 12 4]


#Using LP
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


#complementarity problem
#Not working
# function MCP_Model()
# 	rho_s = [0.2 for i in producers]
# 	kappa_s = [0 0 -7 0]
# 	rho_d = [1 for j in consumers]
# 	kappa_d = [29 46 15 13 42]

# 	m = MCPModel()
# 	@variable(m, x[i in producers, j in consumers] >= 0, start = 100)
# 	@mapping(m, pairs[i in producers, j in consumers], rho_s[i]*sum(x[i, j1] for j1 in consumers) + kappa_s[i] + costs[i,j] + 
# 		rho_d[j]*sum(x[i1, j] for i1 in producers) - kappa_d[j] >= 0)
# 	@complementarity(m, pairs, x)
	
# 	print(m)
# 	solveMCP(m, linear = true)

# 	@show getvalue(x)
# end


LP_Model()
# MCP_Model()