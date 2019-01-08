using JuMP, Complementarity

firms = 1:3
production_process = 1:2
C = [0.55 0.81;
	0.62 1.25;
	0.78 1.35]

K 	= [21000 16000;
	   17000 22000;
	   18000 14000]
#Production capacity
alpha = 2.5 #demand intercept
beta = 0.0000166666667 #inverse demand slope



m = MCPModel()
@variable(m, x[i in firms, b in production_process] >= 0)
@variable(m, q >= 0)
@variable(m, gamma[i in firms, b in production_process] >= 0)#Dual of capacity constraint
@variable(m, lambda) #price
###
@mapping(m, PriceDef, lambda - (alpha - beta*q))
@mapping(m, Mx[i in firms, b in production_process], C[i,b] - lambda + beta*sum(x[i,j] for j in production_process) + gamma[i,b])
@mapping(m, Capacity[i in firms, b in production_process], -x[i,b] + K[i,b])
@mapping(m, MarketClearing, q - sum(x[i,b] for i in firms, b in production_process))

@complementarity(m, PriceDef, q)
@complementarity(m, Mx, x)
@complementarity(m, Capacity, gamma)
@complementarity(m, MarketClearing, lambda)

solution = solveMCP(m, linear = true)

q = getvalue(q); x = getvalue(x)
Price = alpha - beta*q
ConsumerSurplus = alpha*q - 0.5*beta*q^2 - Price*q
@show q
@show Price
@show ConsumerSurplus
