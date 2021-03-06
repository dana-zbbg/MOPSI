# trace amplitude en fonction du pas de temps h

N = 6
Q0 = initial_positions[1:3*N] #sun and first planet
P0 = initial_moments[1:3*N]
N_amp = 20 #nb points pour tracer amplitude
total_time = 10000 #doit suffire pour avoir au moins une periode
list_h = [i*5 for i=1:N_amp]
list_Nh = zeros(Int, N_amp)
for i=1:N_amp
    list_Nh[i] = floor(Int,total_time/(list_h[i]))
end




# ------- Simplectic -------

function plot_amplitude_symp(N_amp)
    list_amplitudes = zeros(N_amp)
    list_log_h = zeros(N_amp)

    for i=1:N_amp
        h = list_h[i]
        Nh = list_Nh[i]
        QN, PN = Symplectic(h, Nh, N, Q0, P0)
        H = energy_H(PN,QN,Nh)
        list_amplitudes[i] = log(maximum(H) - minimum(H))
        list_log_h[i] = log(list_h[i])
    end
    print((list_amplitudes[10]-list_amplitudes[1])/(list_log_h[10]-list_log_h[1]),"\n")

    plot(list_log_h, list_amplitudes, label = "amplitude en fonction du carré du temps", marker =".")
    legend()
    title("Symplectic")
end

#plot_amplitude_symp(N_amp)

# --------- Verlet ---------

function plot_amplitude_verlet(N_amp)
    list_amplitudes = zeros(N_amp)
    list_log_h = zeros(N_amp)

    for i=1:N_amp
        h = list_h[i]
        Nh = list_Nh[i]
        QN, PN = Verlet(h, Nh, N, Q0, P0)
        H = energy_H(PN,QN,Nh)
        list_amplitudes[i] = log(maximum(H) - minimum(H))
        list_log_h[i] = log(list_h[i])
    end
    print((list_amplitudes[10]-list_amplitudes[1])/(list_log_h[10]-list_log_h[1]),"\n")

    plot(list_log_h, list_amplitudes, label = "amplitude en fonction du carré du temps", marker =".")
    legend()
    title("Verlet")
end

#plot_amplitude_verlet(N_amp)


function plot_amplitude_modified_energy(N_amp)
    list_amplitudes = zeros(N_amp)
    list_log_h = zeros(N_amp)

    for i=1:N_amp
        h = list_h[i]
        Nh = list_Nh[i]
        QN, PN = Symplectic(h, Nh, N, Q0, P0)
        H = analytic_EE_modified_energy_H(PN, QN, h, Nh,N)
        list_amplitudes[i] = log(maximum(H) - minimum(H))
        list_log_h[i] = log(list_h[i])
    end
    print((list_amplitudes[10]-list_amplitudes[1])/(list_log_h[10]-list_log_h[1]),"\n")
    
    plot(list_log_h, list_amplitudes, label = "log amplitude modified energy", marker =".")
    xlabel("log time step")
    legend()
    title("Symplectic")
end



function plot_amplitude_modified_energy_V(N_amp)
    list_amplitudes = zeros(N_amp)
    list_log_h = zeros(N_amp)

    for i=1:N_amp
        h = list_h[i]
        Nh = list_Nh[i]
        QN, PN = Verlet(h, Nh, N, Q0, P0)
        H = analytic_modified_energy_H(PN, QN, h, Nh,N)
        list_amplitudes[i] = log(maximum(H) - minimum(H))
        list_log_h[i] = log(list_h[i])
    end
    print((list_amplitudes[10]-list_amplitudes[1])/(list_log_h[10]-list_log_h[1]),"\n")

    plot(list_log_h, list_amplitudes, label = "log amplitude modified energy", marker =".")
    xlabel("log time step")
    legend()
    title("Verlet")
end


