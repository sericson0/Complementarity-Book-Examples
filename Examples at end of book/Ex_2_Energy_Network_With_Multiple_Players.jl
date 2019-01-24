using Complementarity, JuMP, LinearAlgebra
##
#rows for node, columns for producers


function multimarket_trade(nodes, line_data, firm_data, demand_data)
	#get parameters from input dictionaries
	tau = line_data["tau"]; lines_out = line_data["lines_out"]; lines_in = line_data["lines_in"]; line_size = line_data["line_size"]
	C = firm_data["C"]; q_max = firm_data["q_max"]; firm_location = firm_data["firm_location"]
	firms = 1:length(C)
	a = demand_data["a"]; b = demand_data["b"]
	#____________________________________________
	#parameters are
	#		tau_i,j			- marginal cost of transportation between i and j
	#		lines_in 		- dictionary where values are regions which can sell to key
	#		lines_out		- dictionary where values are regions which key can sell to
	#		line_size		- maximum trade from i to j
	#		C 				- marginal costs for each firm
	#		q_max 			- maximum capacity for each firm
	#		firm_location	- the node where each firm is located
	#		a 				- intercept of linear inverse demand curve P = a - b * Q
	#		b 				- slope of linear inverse demand curve P = a - b * Q
	#___________________________________________
	#variables are:
	#		q_f 			- output from each firm
	#		lambda_f		- shadow value of capacity constraint from each firm
	#		delta_i,j 		- trade from node i to node j
	#		domestic_i		- domestic consumption in node i 
	#						(could just allow trade from i to i and combine delta and domestic, but the solver seems to prefer this formulation)
	#		epsilon_i,j 	- shadow value of line constraint when trading from i to j	
	#________________________________________________________________________________________
	#________________________________________________________________________________________
	m = MCPModel()
	@variable(m, q[f in firms] >= 0)
	@variable(m, lambda[f in firms] >= 0)
	@variable(m, delta[i in nodes, j in lines_out[i]] >= 0)
	@variable(m, epsilon[i in nodes, j in lines_out[i]] >= 0)
	@NLexpression(m, domestic[i in nodes], sum(q[f] for f in firms if firm_location[f] == i) - sum(delta[i,j] for j in lines_out[i]))
	#For some reason when I use the expression inside the mappings insdead of hardcoding the price it messes up the solver. Here I save the expression for outputting
	@NLexpression(m, prices[i in nodes], a[i] - b[i]*(domestic[i] + sum(delta[j, i] for j in lines_in[i])))
	#_____________________________________________
	#The complementarity equations say the following:
	#		zero_profit		-This assumes a competitive market and says MC + capacity_value - P >= 0
	#						 Either marginal profit is pushed to zero --after subtracting the shadow value of capacity-- or you will produce zero.
	#		capacity_const  -Either firms produce to their capacity constraint (q == q_max) or the value of additional capacity is zero.
	#		line_constraint -Either the line constraint is met (delta[i,j] = line_size[i,j]) or the value of additional line capacity is zero.
	#		no_arbitrage	-Either you are indifferent to selling domestically and selling to your neighbor (P[i] == P[j] + tau[i,j]) or you don't sell to the neighbor.		
	#_____________________________________________
	@mapping(m, zero_profit[f in firms], C[f] + lambda[f] - (a[firm_location[f]] - b[firm_location[f]]*(domestic[firm_location[f]] + sum(delta[j, firm_location[f]] for j in lines_in[firm_location[f]]))))
	@mapping(m, capacity_const[f in firms], q_max[f] - q[f])
	@mapping(m, line_constraint[i in nodes, j in lines_out[i]], line_size[i][j] - delta[i, j])
	@mapping(m, no_arbitrage[i in nodes, j in lines_out[i]], (a[i] - b[i]*(domestic[i] + sum(delta[j, i] for j in lines_in[i]))) + tau[i,j] + epsilon[i,j] - (a[j] - b[j]*(domestic[j] + sum(delta[j1, j] for j1 in lines_in[j]))))
	##
	#Because we are subsetting delta and epsilon variables to only valid neighbors (not all regions can trade, and trade is not always bidirectional), these values are saved in a JuMP Dict
	#Jump dicts cannot be implicitly iterated through so we must dirctly iterate over these variables, which is why they are in the for loop
	##
	for i in nodes, j in lines_out[i]
		@complementarity(m, line_constraint[i,j], epsilon[i,j])
		@complementarity(m, no_arbitrage[i,j], delta[i,j])
	end
	#
	@complementarity(m, zero_profit, q)
	@complementarity(m, capacity_const, lambda)
	#Lets solve this shit!
	solution = solveMCP(m, linear = true) 
	#Spit out values
	@show getvalue(q)
	@show getvalue(lambda)
	@show getvalue(domestic)
	@show getvalue(delta)
	@show getvalue(prices)
end



nodes = 1:2
N = length(nodes)
line_data = Dict("tau" => 0.5 .*(ones(N) .- Diagonal(ones(N))), "lines_out" => Dict(1 => [2], 2 => []), 
	"lines_in" => Dict(1 => [], 2 => [1]), "line_size" => Dict(1 => [0 5], 2 => [0 0]))
firm_data = Dict( "C" => [10 12 15 18], "q_max" => [10 10 5 5], "firm_location" => [1 1 2 2])
##
#Taking inverse demand of original problem (Didin't realize example used regular demand curve until after coding)
#Original problem is Q[i] = [20 40] - [1 2] P[i]
#Changed to P[i] = [20 40]/[1 2] - Q[i]/[1 2] => P[i] = [20 20] - [1 .5]Q[i]
demand_data = Dict("a" => [20 20], "b" => [1 .5])
##

multimarket_trade(nodes, line_data, firm_data, demand_data)