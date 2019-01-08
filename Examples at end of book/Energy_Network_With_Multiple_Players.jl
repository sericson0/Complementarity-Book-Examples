using Complementarity, JuMP, LinearAlgebra
##
#rows for node, columns for producers
tau_12reg = .5
gamma = [10 12; 
		 15 18]

a = [20; 40]
b = [1; 2]
q_max = [10 10;
		 5 5]
g_12_max = 5
gamma_TSO = 1
nodes = 1:length(a)
prod_per_node = 1:2
##
m = MCPModel()

# println("test1")
@variable(m, s[i in nodes, j in prod_per_node] >= 0)
@variable(m, q[i in nodes, j in prod_per_node] >= 0)
# println("test2")
@variable(m, f[i in prod_per_node] >= 0)
@variable(m, lambda[i in nodes, j in prod_per_node] >= 0)
@variable(m, g_12 >= 0)
@variable(m, eps_12 >= 0)
##
@variable(m, delta[i in nodes, j in prod_per_node])
@variable(m, pi[i in nodes])
# println("test3")
@variable(m, tau_12)
###

@mapping(m, Eq_s[i in nodes, j in prod_per_node], -pi[i] + delta[i,j])
@mapping(m, Eq_q[i in nodes, j in prod_per_node], gamma[i,j] + lambda[i,j] - delta[i,j])
@mapping(m, Eq_f[j in prod_per_node], -pi[1] + (tau_12reg + tau_12) + delta[1,j] + lambda[1,j] - delta[1,j])
# println("test4")
@mapping(m, Eq_lambda[i in nodes, j in prod_per_node], q_max[i,j] - q[i,j])
@mapping(m, Eq_eq[i in nodes, j in prod_per_node], s[i,j] - q[i,j] + ifelse(i == 1, f[j], 0))
# @mapping(m, Eq_eq[i in nodes, j in prod_per_node], s[i,j] - q[i,j] )
# @mapping(m, Eq_eq[i in nodes, j in prod_per_node], s[i,j] - q[i,j] )
@mapping(m, Eq_MC[i in nodes], sum{s[i,j], j in prod_per_node} - a[i] - b[i]*pi[i])
# println("test5")
@mapping(m, Eq_g_12, -tau_12reg - tau_12 + gamma_TSO + eps_12)
@mapping(m, Eq_eps_12, g_12_max - g_12)
@mapping(m, Eq_tau_12, g_12 - sum{f[i], i in prod_per_node})
##
##
@complementarity(m, Eq_s, s)
@complementarity(m, Eq_q, q)
@complementarity(m, Eq_f, f)
@complementarity(m, Eq_lambda, lambda)
@complementarity(m, Eq_g_12, g_12)
@complementarity(m, Eq_eps_12, eps_12)
# println("test7")

print(m)
solution = solveMCP(m, linear = true)

@show getvalue(q)
@show getvalue(lambda)
