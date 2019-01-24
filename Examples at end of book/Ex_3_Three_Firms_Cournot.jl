using JuMP, Complementarity

firms = 1:3
production_process = 1:2
#Matricies defined similar to Matlab, with space indicating columns and ; indicating rows
C = [0.55 0.81;
	0.62 1.25;
	0.78 1.35]

K 	= [21000 16000;
	   17000 22000;
	   18000 14000]
#Production capacity
alpha = 2.5 #demand intercept
beta = 0.0000166666667 #inverse demand slope


#We can see the general setup of MCP models here. 
#	1. Define variables
#	2. For each variable there will be a related complementarity condition
#	3. pair the variables and compelmentarity conditions.
##
#	If you have trouble understanding which variables to pair with wich conditions, it may help to remember that price variables will pair with
#	quantity conditions and quantity variables will pair with price conditions
#	In this example the variables are x (production by firm, process), q (total quantity), gamma (shadow value of capacity), and lambda (price)
#	Hence we have two quantity variables -- x and q -- and two price variables -- gamma and lambda.
#
#	The complementarity constraints are the demand equation P = alpha - beta Q
#	The zero profit condition (marginal costs must be at lest equal to the price once we incorporate the value of capacity)
#	The capacity condition (cant produce more than available capacity)
#	the market clearing condition (must produce as much as is demanded)
# 	Therefore we have two quantity constraints -- capacity constraint and market clearing -- and two price constraints -- demand condition on price and zero profit condition.
#	
#	The pairings are:
#		- demand 	   _|_ q 	  (price constraint with quantity variable)
#		- zero profit  _|_ x 	  (price constraint with quantity variable)
#		- Capacity     _|_ gamma  (quantity constraint with price variable)
#		- Market Clear _|_ lambda (quantity constraint with price variable)
m = MCPModel()
@variable(m, x[i in firms, b in production_process] >= 0)
@variable(m, q >= 0)
@variable(m, gamma[i in firms, b in production_process] >= 0) #Dual of capacity constraint
@variable(m, lambda) #price
###
@mapping(m, DemandEq, lambda - (alpha - beta*q))
@mapping(m, ZeroProfit[i in firms, b in production_process], C[i,b] - lambda + beta*sum(x[i,j] for j in production_process) + gamma[i,b])
@mapping(m, Capacity[i in firms, b in production_process], -x[i,b] + K[i,b])
@mapping(m, MarketClearing, q - sum(x[i,b] for i in firms, b in production_process))

@complementarity(m, DemandEq, q)
@complementarity(m, ZeroProfit, x)
@complementarity(m, Capacity, gamma)
@complementarity(m, MarketClearing, lambda)

print(m)
solution = solveMCP(m, linear = true)

q = getvalue(q); x = getvalue(x)
Price = alpha - beta*q
ConsumerSurplus = alpha*q - 0.5*beta*q^2 - Price*q
@show q
@show Price
@show ConsumerSurplus
