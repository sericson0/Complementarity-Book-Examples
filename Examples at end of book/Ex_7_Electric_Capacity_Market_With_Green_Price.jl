using VariationalInequality, JuMP
#Uses the VariationalInequality package to solve variational inequality problems. 


firms = 1:2
gen_types = ["exist", "convent", "green"]
# type_id = 1:length(gen_types)
#Here I combine all generator data into a dictionary of dictionaries
gen_data = Dict("exist"   => Dict("Mult"=> 1, "MaxCap" =>[140 100], "MC_Const" => [ 180  110], "MCostSlope" => [ 0  0]),
				"convent" => Dict("Mult"=> 1, "MaxCap" =>[ 90  75], "MC_Const" => [1000  900], "MCostSlope" => [15 18]),
				"green"	  => Dict("Mult"=>1.25,"MaxCap"=>[ 50  80], "MC_Const" => [1800 2000], "MCostSlope" => [90 70]))

CapReq = 400

m = VIPModel()
@variable(m, x[h in gen_types, f in firms] >= 0)
@variable(m, mu >= 0)

@constraint(m, CapLim[h in gen_types, f in firms], gen_data[h]["MaxCap"][f] - x[h,f] >= 0)
@mapping(m, CostMult[h in gen_types, f in firms], 
	gen_data[h]["MC_Const"][f] + gen_data[h]["MCostSlope"][f]*x[h, f] - gen_data[h]["Mult"]*mu)

@mapping(m, Capacity, -CapReq + sum(x[h, f] for h in gen_types, f in firms))


for h in gen_types, f in firms
	@innerproduct(m, CostMult[h, f], x[h,f])
end
@innerproduct(m, Capacity, mu)

# print(m)
solveVIP(m, algorithm=:fixed_point, max_iter=10000, tolerance = 1e-3)


Mu =  getvalue(mu)
X  = getvalue(x)
TotPmt = sum(gen_data[h]["Mult"]*Mu*X[h,f] for h in gen_types, f in firms)
TotCostOut = sum(gen_data[h]["MC_Const"][f]*X[h, f] *0.5*gen_data[h]["MCostSlope"][f]*X[h,f]^2 for h in gen_types, f in firms)

CapByType = [sum(X[h, f] for f in firms) for h in gen_types ] 
CostByType = Dict(zip(gen_types, [sum(gen_data[h]["MC_Const"][f] * X[h,f] + 0.5 * gen_data[h]["MCostSlope"][f] * X[h, f]^2 for f in firms) for h in gen_types ]))
ProfitByType = [sum(gen_data[h]["Mult"] * Mu * X[h, f] for f in firms) - CostByType[h] for h in gen_types]

#Compare with table 5.4 on page 214
println("x: ", X, "\nmu: ", Mu, "\nTotPmt: ", TotPmt, "\nTotCostOut: ",TotCostOut, "\nCapByType: ", CapByType, "\nCostByType: ", CostByType, "\nProfitByType: ", ProfitByType)

