# #Code for example problem 1.2.1 -- Three-Variable MCP. I have generalized the problem to 5 players

#The MCP problem here N suppliers competing via a cournot game and a demand curve of alpha - beta*Q.
#There is also a total quantity constraint (say a line constraint) such that sum(q_i) <= Q_max
#The equilibrium is defined by all producers playing optimally (maximizing profits), and total production not exceeding Q_max
#
#The firms solve for their first order condition, where 
#Each firm's first order condition is alpha - beta_i - sum_j{beta_j} - costs_i - lambda_Q
#Here lambda_Q is the shadow value for the production constraint

using Complementarity, JuMP, LinearAlgebra

m = MCPModel()
##
#costs = [1 2] #original problem
costs = [1 1.25 1.5 1.75 2]
beta = 2
#beta = 5
alpha = 10
#Q_max = 1.06
Q_max = 3.25
N = length(costs)
#M is the matrix [2 1 ... 1
# 				  1 2 ... 1
#				  1 1 ... 2]

M = Diagonal(ones(N)) .+ 1
# #
prod_id = 1:N
@variable(m, generation[i in prod_id] >= 0) 
@variable(m, lambda_Q >= 0)
##
@mapping(m, eq[i in prod_id], sum{beta*M[i,j]*generation[i], j in prod_id} - alpha + costs[i] + lambda_Q)
@mapping(m, Q_constraint, Q_max - sum{generation[i], i in prod_id})

@complementarity(m, eq, generation)
@complementarity(m, Q_constraint, lambda_Q)

print(m) 

solution = solveMCP(m, linear = true)

@show getvalue(generation)
@show getvalue(lambda_Q)
