#Code for the Stackelberg MPEC with Firm 2 as leader
using JuMP, Ipopt, Complementarity, Suppressor

function coalModel(C, K, w, Lim, alpha, beta)
	@suppress begin
		firms = 1:size(C,1)
		process = 1:size(C,2)
		m = MCPModel()
		@variable(m, x[i in firms, b in process] >= 0)
		@variable(m, q >= 0)
		@variable(m, gamma[i in firms, b in process] >= 0)
		@variable(m, phi >= 0)
		@variable(m, lamda)
		##
		@mapping(m, PriceDef, lamda - (alpha - beta*q))
		@mapping(m, Mx[i in firms, b in process], C[i,b] - lamda +  beta*sum(x[i,b1] for b1 in process) + gamma[i,b] + w[i]*phi)
		@mapping(m, Capacity[i in firms, b in process], K[i,b] - x[i,b])
		@mapping(m, MarketsClear, q - sum(x[i,b] for i in firms, b in process))
		@mapping(m, YardLimit, Lim - sum(x[i,b] for i in firms, b in process))

		@complementarity(m, PriceDef, q)
		@complementarity(m, Mx, x)
		@complementarity(m, Capacity, gamma)
		@complementarity(m, MarketsClear, lamda)
		@complementarity(m, YardLimit, phi)

		
		solution = solveMCP(m, linear = true)
	X = ones(size(C))
	for i in firms, b in process
		X[i,b] = getvalue(x[i,b])
	end
	return getResults(X, getvalue(q), getvalue(gamma), getvalue(phi), getvalue(lamda), C, alpha, beta)
	end
end



function getResults(x, q, gamma, phi, lamda, C, alpha, beta)
	println(7)
	Price = alpha - beta*q
	ConPmt = Price * q
	ConSurp = alpha*q - 0.5*beta*q^2 - Price*q
	Rev = Price .* sum(x, dims = 2)
	VarCost = sum(C .* x, dims = 2)
	ProdSurp = Rev - VarCost
	TotProdSurp = sum(ProdSurp)
	println(8)
	SocWel = ConSurp + TotProdSurp
	return(Dict("x"=>x, "Rev" => Rev, "VarCost" => VarCost, "ProdSurp" => ProdSurp, "q"=> q, "Price" => Price, "ConPmt" => ConPmt, 
		"ConSurp" => ConSurp, "TotProdSurp" => TotProdSurp, "SocWel" => SocWel, "phi"=>phi))
end

C = [ 0.55 0.81;
	  0.62 1.25;
	  0.78 1.35]

K = [21000 16000;
	 17000 22000;
	 18000 14000]

w = [1 2 2]

Lim = 60000
alpha = 2.5
beta = 0.0000166666667

results = coalModel(C, K, w, Lim, alpha, beta)

# println(results)
for (key, val) in results
	println(key, ": ", val)
end