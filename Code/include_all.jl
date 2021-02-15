# INITIAL CONDITIONS
include("variables_initial_conditions")

# GRADIENT AND F
include("potential_V_function_f.jl")

# NUMERIC SCHEMES
include("Explicite_Euler.jl")
include("Implicite_Euler.jl")
include("RK_4.jl")
include("Symplectic.jl")
include("Verlet.jl")
include("Verlet_B.jl")

# ENERGY
include("energy_functions.jl")
include("amplitude.jl")

# DISPLAY
include("display_functions.jl")
include("animated_plot.jl")
include("energy_plot.jl")
