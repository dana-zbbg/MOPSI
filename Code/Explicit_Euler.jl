##Explicit Euler Scheme
##solve by explicit euler for time span h, Nh iterations,
##initialised at position Q0 at movement quantity P0
##Q0[3*i-2:3*i] = initial position of particule i (idem P0)
function Explicit_Euler(h, Nh, N, Q0, P0)
    QN = zeros(3*N,Nh) #positions: Q[3*i-2:3*i,j] = composantes de la position de la particule i à iteration j
    PN = zeros(3*N,Nh) #quantite de mouvement tq P[3*i-2:3*i,j] = composantes de quantite mvt de la particule i à interation j
    Y = zeros(6*N) #vector (Q,P) at time n
    Y[1:3*N] = Q0
    Y[3*N+1:6*N] = P0
    QN[1:3*N,1] = Q0
    PN[1:3*N,1] = P0
    for i = 2:Nh
        Y = Y+h*f(Y, N) ##euler explicite
        QN[:,i] = Y[1:3*N] ##storage of positons at time i in QN
        PN[:,i] = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    return QN, PN
end
