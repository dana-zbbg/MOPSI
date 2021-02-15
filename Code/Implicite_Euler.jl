##Euler Implicite Scheme
##solve by explicit euler for time span h, Nh iterations,
##initialised at position Q0 at movement quantity P0
function Implicite_Euler(h, Nh, N, Q0, P0, eps = 10e-4)
    QN = zeros(3*N,Nh) #positions tq Q[3*i-2:3*i,j] = composantes de la position de la particule i à iteration j
    PN = zeros(3*N,Nh) #quantite de mouvement tq P[3*i-2:3*i,j] = composantes de quantite mvt de la particule i à interation j
    Y = zeros(6*N) #vector (Q,P) at time n
    Y[1:3*N] = Q0
    Y[3*N+1:6*N] = P0
    Z0 = zeros(6*N)
    Z1 = zeros(6*N)
    for i= 1:Nh
        ## we find Z1 : an approximation of Y at time i+1
        println(i)
        ZInit = Y + h*f(Y, N)
        Stop_condition = eps*norm(ZInit)
        Z0 = ZInit
        Z1 = Y+h*f(Z0, N)
        j=0
        while (norm(Z1-Z0)> Stop_condition && j<100)
            Z0 = Z1
            Z1 = Y+h*f(Z0, N)
            j+=1
            print("iter ")
            println(j)
        end
        Y = Y+h*f(Z1, N) ##euler implicite
        QN[:,i] = Y[1:3*N] ##storage of positons at time i in QN
        PN[:,i] = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    return QN, PN
end
