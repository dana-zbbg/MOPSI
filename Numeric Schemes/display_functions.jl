function display(QN, N)
    figure()
    for i=0:N-1
        x = QN[1+i*3,:]
        y = QN[2+i*3,:]
        z = QN[3+i*3,:]
        plot3D(x,y,z, label = particules[i+1], color = colors[i+1], linewidth=0.8)
    end
    PyPlot.legend()
    #scatter3D(x,y,z, c=z)
end
