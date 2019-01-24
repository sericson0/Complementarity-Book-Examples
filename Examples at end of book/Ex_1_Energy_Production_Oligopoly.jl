# #Code for example problem 1.2.5 -- Energy Duopoly MCP. I have also generalized the problem to 5 players

#The MCP problem here N suppliers competing via a cournot game and a demand curve of alpha - beta*Q.
#There is also a total quantity constraint (say a transmission line constraint) such that sum(q_i) <= Q_max
#The equilibrium is defined by all producers playing optimally (maximizing profits), and total production not exceeding Q_max
#
#The firms solve for their first order condition, where 
#Each firm's first order condition is alpha - beta_i - sum_j{beta_j} - costs_i - lambda_Q
#Here lambda_Q is the shadow value for the production constraint

#Using loads packages, and is equivalent to library() in R or import in Python
#You will need to install these packages before running the code. 
using Complementarity, JuMP, LinearAlgebra

#Creates a user defined function to solve an oligopoly equilibrium using an MCP framework
function oligopolyMCP(costs, alpha, beta, Q_max, print_equation = true) 
	#see https://github.com/chkwon/Complementarity.jl for description of MCP modeling in Julia
	m = MCPModel()
	N = length(costs)
	M = Diagonal(ones(N)) .+ 1
	#M is the matrix [2 1 ... 1
	# 				  1 2 ... 1
	#				  1 1 ... 2]
	prod_id = 1:N
	#we first define the variables in the model and can supply upper and lower bounds
	@variable(m, generation[i in prod_id] >= 0) 
	@variable(m, lambda_Q >= 0)
	#@mapping determines the equations for the complementarity problems. The second input (eq and Q_constraint)
	#names the complementarity equation while the third input defines the complementarity equation 
	@mapping(m, eq[i in prod_id], sum(beta*M[i,j]*generation[j] for j in prod_id) - alpha + costs[i] + lambda_Q)
	@mapping(m, Q_constraint, Q_max - sum{generation[i], i in prod_id})
	#@complementarity connects equations to variables. 
	@complementarity(m, eq, generation)
	@complementarity(m, Q_constraint, lambda_Q)
	#Prints model formulation
	if print_equation == true
		print(m) 
	end
	#solve model
	solution = solveMCP(m, linear = true)
	#display results
	@show getvalue(generation)
	@show getvalue(lambda_Q)
end
###

#Book example
costs = [1 2] 
alpha = 10
beta = 5
Q_max = 1.06
oligopolyMCP(costs, alpha, beta, Q_max)

#Increase to 5 firms
costs = [1 1.25 1.5 1.75 2]
alpha = 10
beta = 2
Q_max = 3.25
oligopolyMCP(costs, alpha, beta, Q_max)

