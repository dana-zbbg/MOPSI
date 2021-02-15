include("include_all.jl")

h=10
Nh=4000

time_array = [h*i for i=1:Nh]

Q0 = initial_positions[1:3*N] #sun and first planet
P0 = initial_moments[1:3*N]


# Verlet
QN_Verlet, PN_Verlet = Verlet(h, Nh, N, Q0, P0)
H_Verlet = energy_H(PN_Verlet, QN_Verlet, Nh)
plot(time_array, H_Verlet, label = "Verlet")
H_Verlet_modified = analytic_modified_energy_H(PN_Verlet, QN_Verlet,h,Nh)
plot(time_array, H_Verlet_modified, label = "Verlet modified energy")

#=
# Symplectic
QN_Sym, PN_Sym = Symplectic(h, Nh, N, Q0, P0)
H_Sym = energy_H(PN_Sym, QN_Sym, Nh)
plot(time_array, H_Sym, label = "Symplectic")
H_Symplectic_modified = analytic_modified_energy_H(PN_Verlet, QN_Verlet,h,Nh)
plot(time_array, H_Symplectic_modified, label = "Symplectic modified energy")



# Implicite
QN_Imp, PN_Imp = Implicite_Euler(h, Nh, N, Q0, P0)
H_Imp = energy_H(PN_Imp, QN_Imp, Nh)
plot(time_array, H_Imp, label = "Implicit Euler")


# Implicite_milieu
QN_Imp_m, PN_Imp_m = implicite_milieu(h, Nh, N, Q0, P0)
H_Imp_m = energy_H(PN_Imp_m, QN_Imp_m, Nh)
plot(time_array, H_Imp_m, label = "Implicite milieu")


# Explicit
QN_Exp, PN_Exp = Explicite_Euler(h, Nh, N, Q0, P0)
H_Exp = energy_H(PN_Exp, QN_Exp, Nh)
plot(time_array, H_Exp, label = "Explicit Euler")


# RK4
QN_RK, PN_RK = RK_4(h, Nh, N, Q0, P0)
H_RK = energy_H(PN_RK, QN_RK, Nh)
plot(time_array, H_RK, label = "Runge-Kutta 4")
=#
xlabel("time (in days)")
ylabel("Energy")
legend()
