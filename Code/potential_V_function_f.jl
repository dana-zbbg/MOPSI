## gradient_V computes the gradient of V at point q_k evaluated for position vector Q, masses M
function gradient_V(k, Q, N)
    S = zeros(3)
    vect = zeros(3)
    for i = 1:N
        if (i != k)
            vect = Q[3*k-2:3*k]-Q[3*i-2:3*i]
            S += masses[k]*masses[i]*(1/norm(vect)^3)*vect
        end
    end
    return S*G
end


##\gradient_q H
function gra_q_H(Q)
    K=N
    Gra_q=zeros(size(Q))
    for k=1:K
        Gra_q[(3*k-2):(3*k)]=gradient_V(k, Q, N)
    end
    return Gra_q
end


##\gradient_p H
function gra_p_H(P)
    K=N
    Gra_p=zeros(size(P))
    for k=1:K
        Gra_p[(3*k-2):(3*k)]=1/masses[k]*P[(3*k-2):(3*k)]
    end
    return Gra_p
end


## function of the differential system :
##Y' = f(Y, N, M, MN) where Y = (Q,P)
function f(Y, N)
    Q = Y[1:3*N]
    Z = zeros(6*N)
    for k = 1:N #planete 1 to N
        Z[3*k-2 : 3*k] = Y[3*(N+k)-2 : 3*(N+k)]/masses[k]
        Z[3*(N+k)-2 : 3*(N+k)] = -gradient_V(k, Q, N)[:]
    end
    return Z
end
