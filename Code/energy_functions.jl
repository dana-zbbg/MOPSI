#energy_H computes the total energy of the system
function energy_H(P,Q,Nh,N)
    H=zeros(Nh)
    for k=1:Nh
        for i=1:N
            H[k]=H[k]+0.5*(1/masses[i])*P[3*i-2:3*i,k]'*P[3*i-2:3*i,k]
        end
        for i=2:N
            for j=1:(i-1)
                H[k]=H[k]-G*masses[i]*masses[j]/norm(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])
            end
        end
    end
    return H

end

## it is just a part of calculation of the approximation of modified hamitonien
function U_energy(Q,k,N)
    E=0;
    for i=2:N
        for j=1:(i-1)
            E=E-G*masses[i]*masses[j]/norm(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])
        end
    end
    return E
end

#modified_energy_H
function modified_energy_H(P,Q,Nh,N)
    H=energy_H(P,Q,Nh,N)
    Qneg=zeros(3*N,1)
    Pneg=zeros(3*N,1)
    Y=zeros(6*N,1)
    Y[1:3*N] = Q[:,1]
    Y[3*N+1:6*N] = P[:,1]
    h1=-h
    for i= 1:Nh
        Y[3*N+1:6*N]=Y[3*N+1:6*N]-h1/2*gra_q_H(Y[1:3*N],N)
        Y[1:3*N]=Y[1:3*N]+h1*gra_p_H(Y[3*N+1:6*N],N)
        Y[3*N+1:6*N]=Y[3*N+1:6*N]-h1/2*gra_q_H(Y[1:3*N],N)

        Qneg = Y[1:3*N] ##storage of positons at time i in QN
        Pneg = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    k=1
    Temp1=-gra_q_H(Q[:,k+1],N)
    Temp2=-gra_q_H(Qneg,N)
    Temp3=-gra_q_H(Q[:,k],N)
    H[1]=H[1]+1/4*(U_energy(Q,2,N)-2*U_energy(Q,1,N)+U_energy(Qneg,1,N))
    for i=1:N
            if i!=2
                H[k]=H[k]+1/6*h*P[3*i-2:3*i,k]'*1/masses[i]*0.5*(Temp1[3*i-2:3*i]-Temp2[3*i-2:3*i])
                H[k]=H[k]+5/24*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*Temp3[3*i-2:3*i]
                H[k]=H[k]+1/12*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*(Temp1[3*i-2:3*i]-2*Temp3[3*i-2:3*i]+Temp2[3*i-2:3*i])
            end
    end

    for k=2:Nh-1
        Temp1=-gra_q_H(Q[:,k+1],N)
        Temp2=-gra_q_H(Q[:,k-1],N)
        Temp3=-gra_q_H(Q[:,k],N)
        H[k]=H[k]+1/4*(U_energy(Q,k+1,N)-2*U_energy(Q,k,N)+U_energy(Q,k-1,N))
        for i=1:N
            if i!=2
                H[k]=H[k]+1/6*h*P[3*i-2:3*i,k]'*1/masses[i]*0.5*(Temp1[3*i-2:3*i]-Temp2[3*i-2:3*i])
                H[k]=H[k]+5/24*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*Temp3[3*i-2:3*i]
                H[k]=H[k]+1/12*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*(Temp1[3*i-2:3*i]-2*Temp3[3*i-2:3*i]+Temp2[3*i-2:3*i])
            end
        end
    end
    return H
end


#modified energy for Verlet scheme A
function analytic_modified_energy_H(P,Q,h,Nh,N)
    H=energy_H(P,Q,Nh,N)
    Mneg=zeros(3*N,3*N)
    M=zeros(3*N,3*N)
    V=zeros(3*N)
    for k=1:N
       V[3*k-2:3*k].=masses[k]
    end
    M=Diagonal(V)
    Mneg=inv(M)

    for k=1:Nh
        Temp1=-gra_q_H(Q[:,k],N)
        Heiss=zeros(3*N,3*N)
        I=Diagonal([1,1,1])
        for i=1:N
            for j=1:N
                if i==j
                    for k_1=1:N
                        if k_1 != i
                            dis=norm(Q[3*i-2:3*i,k]-Q[3*k_1-2:3*k_1,k])
                            A=(Q[3*i-2:3*i,k]-Q[3*k_1-2:3*k_1,k])*(Q[3*i-2:3*i,k]-Q[3*k_1-2:3*k_1,k])'
                            Heiss[3*i-2:3*i,3*i-2:3*i]=Heiss[3*i-2:3*i,3*i-2:3*i].+masses[i]*masses[k_1]/dis^3*I
                            Heiss[3*i-2:3*i,3*i-2:3*i]=Heiss[3*i-2:3*i,3*i-2:3*i].-3*masses[i]*masses[k_1]/dis^5*A
                        end
                    end

                else
                    dis=norm(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])
                    A=(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])*(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])'
                    Heiss[3*i-2:3*i,3*j-2:3*j]=Heiss[3*i-2:3*i,3*j-2:3*j].-masses[i]*masses[j]/dis^3*I
                    Heiss[3*i-2:3*i,3*j-2:3*j]=Heiss[3*i-2:3*i,3*j-2:3*j].+masses[i]*masses[j]*3/dis^5*A
                end
            end
        end
        Heiss=G*Heiss

        H[k]=H[k]+h^2*((1/12)*P[:,k]'*(Mneg^1)*Heiss*(Mneg^1)*P[:,k]-(1/24)*Temp1'*Mneg^1*Temp1)
        #print("k: ",k,"  ",P[:,k]'*(Mneg^0)*Heiss*(Mneg^0)*P[:,k],"  ",Temp1'*Mneg^0*Temp1,"\n")
    end
    return H
end


##analytic modified Hamiltonien of symplectic
function analytic_EE_modified_energy_H(P,Q,h,Nh,N)
    H=energy_H(P,Q,Nh,N)
    Mneg=zeros(3*N,3*N)
    M=zeros(3*N,3*N)
    V=zeros(3*N)
    for k=1:N
       V[3*k-2:3*k].=masses[k]
    end
    M=Diagonal(V)
    Mneg=inv(M)

    for k=1:Nh
        Temp1=-gra_q_H(Q[:,k],N)
        H[k]=H[k]+h*(1/2*P[:,k]'*(Mneg^1)*Temp1)
    end
    return H
end



##The main fonction of the approximation of modified hamitonien
function modified_energy_H(P,Q,h,Nh,N)
    H=energy_H(P,Q,Nh,N)
    Qneg=zeros(3*N,1)
    Pneg=zeros(3*N,1)
    Y=zeros(6*N,1)
    Y[1:3*N] = Q[:,1]
    Y[3*N+1:6*N] = P[:,1]
    for i= 1:Nh
        Y[3*N+1:6*N]=Y[3*N+1:6*N]+h/2*gra_q_H(Y[1:3*N],N)
        Y[1:3*N]=Y[1:3*N]-h*gra_p_H(Y[3*N+1:6*N],N)
        Y[3*N+1:6*N]=Y[3*N+1:6*N]+h/2*gra_q_H(Y[1:3*N],N)

        Qneg = Y[1:3*N] ##storage of positons at time i in QN
        Pneg = Y[3*N+1:6*N] ##storage of quantity of movement at time i in PN
    end
    k=1
    Temp1=-gra_q_H(Q[:,k+1],N)
    Temp2=-gra_q_H(Qneg,N)
    Temp3=-gra_q_H(Q[:,k],N)
    H[1]=H[1]+1/4*(U_energy(Q,2,N)-2*U_energy(Q,1,N)+U_energy(Qneg,1,N))
    for i=1:N

            H[k]=H[k]+1/6*h*P[3*i-2:3*i,k]'*1/masses[i]*0.5*(Temp1[3*i-2:3*i]-Temp2[3*i-2:3*i])
            H[k]=H[k]+5/24*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*Temp3[3*i-2:3*i]
            H[k]=H[k]+1/12*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*(Temp1[3*i-2:3*i]-2*Temp3[3*i-2:3*i]+Temp2[3*i-2:3*i])

    end

    for k=2:Nh-1
        Temp1=-gra_q_H(Q[:,k+1],N)
        Temp2=-gra_q_H(Q[:,k-1],N)
        Temp3=-gra_q_H(Q[:,k],N)
        H[k]=H[k]+1/4*(U_energy(Q,k+1,N)-2*U_energy(Q,k,N)+U_energy(Q,k-1,N))
        for i=1:N

            H[k]=H[k]+1/6*h*P[3*i-2:3*i,k]'*1/masses[i]*0.5*(Temp1[3*i-2:3*i]-Temp2[3*i-2:3*i])
            H[k]=H[k]+5/24*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*Temp3[3*i-2:3*i]
            H[k]=H[k]+1/12*h^2*Temp3[3*i-2:3*i]'*1/masses[i]*(Temp1[3*i-2:3*i]-2*Temp3[3*i-2:3*i]+Temp2[3*i-2:3*i])

        end

    end
    return H
end


##modified Hamiltonien of Verlet(position first)
function analytic_modified_energy_H_B(P,Q,Nh,N)
    Mneg=zeros(3*N,3*N)
    M=zeros(3*N,3*N)
    V=zeros(3*N)
    for k=1:N
       V[3*k-2:3*k].=masses[k]
    end
    M=Diagonal(V)
    Mneg=inv(M)
    H=energy_H(P,Q,Nh,N)
    for k=1:Nh
        Temp1=-gra_q_H(Q[:,k],N)
        Heiss=zeros(3*N,3*N)
        I=Diagonal([1,1,1])
        for i=1:N
            for j=1:N
                if i==j
                    for k_1=1:N
                        if k_1 != i
                            dis=norm(Q[3*i-2:3*i,k]-Q[3*k_1-2:3*k_1,k])
                            A=(Q[3*i-2:3*i,k]-Q[3*k_1-2:3*k_1,k])*(Q[3*i-2:3*i,k]-Q[3*k_1-2:3*k_1,k])'
                            Heiss[3*i-2:3*i,3*i-2:3*i]=Heiss[3*i-2:3*i,3*i-2:3*i].+masses[i]*masses[k_1]/dis^3*I
                            Heiss[3*i-2:3*i,3*i-2:3*i]=Heiss[3*i-2:3*i,3*i-2:3*i].-3*masses[i]*masses[k_1]/dis^5*A
                        end
                    end

                else
                    dis=norm(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])
                    A=(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])*(Q[3*i-2:3*i,k]-Q[3*j-2:3*j,k])'
                    Heiss[3*i-2:3*i,3*j-2:3*j]=Heiss[3*i-2:3*i,3*j-2:3*j].-masses[i]*masses[j]/dis^3*I
                    Heiss[3*i-2:3*i,3*j-2:3*j]=Heiss[3*i-2:3*i,3*j-2:3*j].+masses[i]*masses[j]*3/dis^5*A
                end
            end
        end
        Heiss=G*Heiss

        H[k]=H[k]+h^2*((-1/24)*P[:,k]'*(Mneg^1)*Heiss*(Mneg^1)*P[:,k]+(1/6)*Temp1'*Mneg^1*Temp1)
        #print("k: ",k,"  ",P[:,k]'*(Mneg^0)*Heiss*(Mneg^0)*P[:,k],"  ",Temp1'*Mneg^0*Temp1,"\n")

    end

    return H

end
