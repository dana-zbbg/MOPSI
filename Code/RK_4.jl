#RK_4
function RK_4(h, Nh, N, Q0, P0)
    QN = zeros(3*N,Nh) #positions tq Q[3*i-2:3*i,j] = composantes de la position de la particule i à iteration j
    PN = zeros(3*N,Nh) #quantite de mouvement tq P[3*i-2:3*i,j] = composantes de quantite mvt de la particule i à interation j
    Y = zeros(6*N) #vector (Q,P) at time n
    Y[1:3*N] = Q0
    Y[3*N+1:6*N] = P0
    for i_1=1:Nh
        if i_1==1
            Q=Q0
            P=P0
        else
            Q=QN[:,i_1-1]
            P=PN[:,i_1-1]
        end
        k_1=[gra_p_H(P,N);-gra_q_H(Q,N)]
        k_2=[gra_p_H(P+1/3*h*k_1[3*N+1:6*N],N);-gra_q_H(Q+1/3*h*k_1[1:3*N],N)]
        k_3=[gra_p_H(P-1/3*h*k_1[3*N+1:6*N]+h*k_2[3*N+1:6*N],N);-gra_q_H(Q-1/3*h*k_1[1:3*N]+h*k_2[1:3*N],N)]
        k_4=[gra_p_H(P+h*k_1[3*N+1:6*N]-h*k_2[3*N+1:6*N]+h*k_3[3*N+1:6*N],N);-gra_q_H(Q+h*k_1[1:3*N]-h*k_2[1:3*N]+h*k_3[1:3*N],N)]
        Y=Y+h/8*(k_1+3*k_2+3*k_3+k_4)

        QN[:,i_1] = Y[1:3*N] ##storage of positons at time i in QN
        PN[:,i_1] = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    return QN, PN
end
