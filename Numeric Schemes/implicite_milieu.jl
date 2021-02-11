##Euler midpoint Implicite Scheme
function implicte_milieu(h, Nh, N, Q0, P0)
    QN = zeros(3*N,Nh) #positions tq Q[3*i-2:3*i,j] = composantes de la position de la particule i à iteration j
    PN = zeros(3*N,Nh) #quantite de mouvement tq P[3*i-2:3*i,j] = composantes de quantite mvt de la particule i à interation j
    Y = zeros(6*N) #vector (Q,P) at time n
    Y[1:3*N] = P0
    Y[3*N+1:6*N] = Q0
    for i_1=1:Nh
        if i_1==1
            Q=Q0
            P=P0
        else
            Q=QN[:,i_1-1]
            P=PN[:,i_1-1]
        end
        for j=1:7
            Z1=[-gra_q_H(Q);gra_p_H(P)]
            Z2=[-gra_q_H(Y[3*N+1:6*N]);gra_p_H(Y[1:3*N])]

            Y=Y+h/2*(Z1+Z2)
        end

        PN[:,i_1] = Y[1:3*N] ##storage of positons at time i in QN
        QN[:,i_1] = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    return QN, PN
end
