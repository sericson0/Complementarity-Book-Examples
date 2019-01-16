#In this exampleis from 3.4.2.1 and 3.4.2.2
#We solve a welfare maximization problem both as an NLP and as an MCP
#There are three firms, each with two production technologies. Costs C and capacities K are given below
using JuMP, Ipopt, Complementarity, Suppressor

firms = 1:3
production_process = 1:2
C = [0.55 0.81;
	0.62 1.25;
	0.78 1.35]

K = [21000 16000;
     17000 22000;
     18000 14000]
#Production capacity
alpha = 2.5 #demand intercept
beta = 0.0000166666667 #inverse demand slope

#We first solve as an NLP using the IpoptSolver
m = Model(solver = IpoptSolver())

@variable(m, x[i in firms, b in production_process] >= 0)
@variable(m, q >= 0)
@variable(m, NegSocWel)
##
#For LP, NLP, MPEC, ... we define an objective function to Max or Min. Here we minimize negative social welfare
@objective(m, Min, NegSocWel)
#Constraints can be equality or inequality constraints
#Defines NegSocWel
@constraint(m, sum(C[i, b]*x[i,b] for i in firms, b in production_process) - (alpha*q - 0.5*beta*q^2) - NegSocWel == 0)
#Capacity Constraint
@constraint(m, Capacity[i in firms, b in production_process], K[i, b] - x[i, b] >= 0)
#Market clearing condition
@constraint(m, MarketClearing, q - sum(x[i, b] for i in firms, b in production_process) == 0)

@suppress begin
	solve(m)
end
NLP_q = getvalue(q)
#############################################
#############################################
#Now solve as an MCP
m = MCPModel()
#Variables are production (x[i,b]), total quantity (q),
#along with the dual variables for value of capacity for each firm-technology (gamma[i,b]) and the price (lambda)
@variable(m, x[i in firms, b in production_process] >= 0)
@variable(m, q >= 0)
@variable(m, gamma[i in firms, b in production_process] >= 0) #Dual of capacity constraint
@variable(m, lambda) #price
##
#Prices are determined by demand curve or are zero
@mapping(m, PriceDef, lambda - (alpha - beta*q))
#Optimal production decision says marginal costs plus the shadow value of capacity are equal to the price, or production is zero
@mapping(m, Mx[i in firms, b in production_process], C[i,b] - lambda + gamma[i,b])
#Firms produce at full capacity or the value of capacity is zero
@mapping(m, Capacity[i in firms, b in production_process], K[i,b]-x[i,b] )
#The market clears or the price is zero
@mapping(m, MarketClearing, q - sum(x[i,b] for i in firms, b in production_process))

@complementarity(m, PriceDef, q)
@complementarity(m, Mx, x)
@complementarity(m, Capacity, gamma)
@complementarity(m, MarketClearing, lambda)

@suppress begin
	solution = solveMCP(m, linear = true)
end

MCP_q = getvalue(q)
#If all goes according to plan then the two quantities should be the same
println("quantity from NLP: ", NLP_q,", quantity from MCP: ",MCP_q)
