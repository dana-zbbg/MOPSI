##verlet Scheme (position first)
function Verlet_B(h, Nh, N, Q0, P0)
    QN = zeros(3*N,Nh) #positions tq Q[3*i-2:3*i,j] = composantes de la position de la particule i à iteration j
    PN = zeros(3*N,Nh) #quantite de mouvement tq P[3*i-2:3*i,j] = composantes de quantite mvt de la particule i à interation j
    V=zeros(3*N)
    for k=1:N
       V[3*k-2:3*k].=masses[k]
    end
    M=Diagonal(V)
    Mneg=inv(M)
    QN[:,1]=Q0+h/2*Mneg*P0
    PN[:,1]=P0-h*gra_q_H(QN[:,1],N)
    QN[:,1]=QN[:,1]+h/2*Mneg*PN[:,1]
    for i= 2:Nh
        QN[:,i]=QN[:,i-1]+h/2*Mneg*PN[:,i-1]
        PN[:,i]=PN[:,i-1]-h*gra_q_H(QN[:,i],N)
        QN[:,i]=QN[:,i]+h/2*Mneg*PN[:,i]
    end
    return QN, PN
end
