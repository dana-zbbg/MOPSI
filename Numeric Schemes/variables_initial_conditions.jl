using PyPlot
using LinearAlgebra


## ---------- VARIABLES ----------
h=10 #h = time step
Nh=2000 #Nh = number of iterations
N=6 #N = nb of planets

l = 1 #u.a = unite astronomique
m0 = 1 #solar mass

G = 2.95912208286*10^(-4) #constante gravitationnelle

# ------------- DATA --------------
particules = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
colors = ["red","blue","orange","firebrick","green","purple"]
masses = [1.00000597682, 0.000954786104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 1/(1.3*10^8)]
initial_positions = [0., 0., 0., -3.5023653, -3.8169847, -1.5507963, 9.0755314, -3.0458353, -1.6483708, 8.3101420, -16.2901086, -7.2521278, 11.4707666, -25.7294829, -10.8169456, -15.5387357, -25.2225594, -3.1902382]
initial_velocities = [0., 0., 0., 0.00565429, -0.00412490, -0.00190589, 0.00168318, 0.00483525, 0.00192462, 0.00354178, 0.00137102, 0.00055029, 0.00288930, 0.00114527, 0.00039677, 0.00276725, -0.00170702, -0.00136504]
initial_moments = zeros(3*6)
for i=1:6
    initial_moments[3*i-2:3*i] = initial_velocities[3*i-2:3*i]*masses[i]
end
