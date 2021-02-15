include("include_all.jl")
using CSV, DataFrames


h = 20
N = 6
Nh = 5000
Q0 = initial_positions[1:3*N] #sun and first planet
P0 = initial_moments[1:3*N]

#QN = (3*N,Nh) Q[3*i-2:3*i,j] = position particule i à iteration j
QN_Exp, PN_Exp = Explicit_Euler(h, Nh, N, Q0, P0)
#QN_Verlet, PN_Verlet = Verlet(h, Nh, N, Q0, P0)

CSV.write("explicit.csv", DataFrame(QN_Exp), header = false)
CSV.write("verlet.csv", DataFrame(QN_Verlet), header = false)
