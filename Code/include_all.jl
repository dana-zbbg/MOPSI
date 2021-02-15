# INITIAL CONDITIONS
include("variables_initial_conditions.jl")

# GRADIENT AND F
include("potential_V_function_f.jl")

# NUMERIC SCHEMES
include("Explicit_Euler.jl")
include("Implicit_Euler.jl")
include("RK_4.jl")
include("Symplectic.jl")
include("Verlet.jl")
include("Verlet_B.jl")

# ENERGY
include("energy_functions.jl")

# DISPLAY
include("display_functions.jl")
