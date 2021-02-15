include("include_all.jl")

h = 20
Nh = 10000


# ----------- explicit ----------
Q, P = Explicit_Euler(h, Nh, N, Q0, P0)
display(Q,N)
title("Explicit Euler")

# ----------- implicit ----------
Q, P = Implicit_Euler(h, 208, N, Q0, P0) #Nh reduit sinon ca diverge et on ne voit plus rien
display(Q,N)
title("Implicit Euler")

# -------- Runge Kutta 4 --------
Q, P = RK_4(h, Nh, N, Q0, P0)
display(Q,N)
title("Runge Kutta 4")

# ---------- symplectic ---------
Q, P = Symplectic(h, Nh, N, Q0, P0)
display(Q,N)
title("Symplectic Euler")

# ----------- Verlet ------------
Q, P = Verlet(h, Nh, N, Q0, P0)
display(Q,N)
title("Verlet")

# ---------- Verlet B -----------
Q, P = Verlet_B(h, Nh, N, Q0, P0)
display(Q,N)
title("Verlet bis")
