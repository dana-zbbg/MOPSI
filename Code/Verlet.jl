##verlet Scheme
##solve by verlet for time span h, Nh iterations,
##initialised at position Q0 at movement quantity P0
function Verlet(h, Nh, N, Q0, P0)
    QN = zeros(3*N,Nh) #positions tq Q[3*i-2:3*i,j] = composantes de la position de la particule i à iteration j
    PN = zeros(3*N,Nh) #quantite de mouvement tq P[3*i-2:3*i,j] = composantes de quantite mvt de la particule i à interation j
    Y = zeros(6*N) #vector (Q,P) at time n
    Y[1:3*N] = Q0
    Y[3*N+1:6*N] = P0
    for i= 1:Nh
        Y[3*N+1:6*N]=Y[3*N+1:6*N]-h/2*gra_q_H(Y[1:3*N],N)
        Y[1:3*N]=Y[1:3*N]+h*gra_p_H(Y[3*N+1:6*N],N)
        Y[3*N+1:6*N]=Y[3*N+1:6*N]-h/2*gra_q_H(Y[1:3*N],N)

        QN[:,i] = Y[1:3*N] ##storage of positons at time i in QN
        PN[:,i] = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    return QN, PN
end

