import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.mplot3d.axes3d as p3
import os
import csv

path_animation = 'C:/Users/Audre/Documents/Ponts/2A/S1/MOPSI/clean_code/MOPSI/Figures/animated'
path_CSV = 'C:/Users/Audre/Documents/Ponts/2A/S1/MOPSI/clean_code/MOPSI/Ressources'


particules = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
colors = ["red","blue","orange","firebrick","green","purple"]
N = 6
Nh = 5000

N_anim = 1000


file_CSV = open('C:/Users/Audre/Documents/Ponts/2A/S1/MOPSI/clean_code/MOPSI/Ressources/explicit.csv')
data_CSV = csv.reader(file_CSV)
list_CSV = list(data_CSV)


x_maxi = 0
x_mini = 0
y_maxi = 0
y_mini = 0
z_maxi = 0
z_mini = 0
QN = np.zeros((N*3,Nh))
for i in range(N):
    for j in range(Nh):
        QN[3*i,j] = float(list_CSV[3*i][j])
        QN[3*i+1,j] = float(list_CSV[3*i+1][j])
        QN[3*i+2,j] = float(list_CSV[3*i+2][j])
        if (QN[3*i,j]>x_maxi):
            x_maxi = QN[3*i,j]
        elif (QN[3*i,j]<x_mini):
            x_mini = QN[3*i,j]
        if (QN[3*i+1,j]>y_maxi):
            y_maxi = QN[3*i+1,j]
        elif (QN[3*i+1,j]<y_mini):
            y_mini = QN[3*i+1,j]
        if (QN[3*i+2,j]>z_maxi):
            z_maxi = QN[3*i+2,j]
        elif (QN[3*i+2,j]<z_mini):
            z_mini = QN[3*i+2,j]



fig = plt.figure()
ax = p3.Axes3D(fig)



ax.set_xlim3d([x_mini, x_maxi])
ax.set_xlabel('X')

ax.set_ylim3d([y_mini, y_maxi])
ax.set_ylabel('Y')

ax.set_zlim3d([z_mini, z_maxi])
ax.set_zlabel('Z')

ln1, = ax.plot(QN[0,0:1],   QN[1,0:1],     QN[2,0:1],     lw=1, color=colors[0])
ln2, = ax.plot(QN[3*1,0:1], QN[3*1+1,0:1], QN[3*1+2,0:1], lw=1, color=colors[1])
ln3, = ax.plot(QN[3*2,0:1], QN[3*2+1,0:1], QN[3*2+2,0:1], lw=1, color=colors[2])
ln4, = ax.plot(QN[3*3,0:1], QN[3*3+1,0:1], QN[3*3+2,0:1], lw=1, color=colors[3])
ln5, = ax.plot(QN[3*4,0:1], QN[3*4+1,0:1], QN[3*4+2,0:1], lw=1, color=colors[4])
ln6, = ax.plot(QN[3*5,0:1], QN[3*5+1,0:1], QN[3*5+2,0:1], lw=1, color=colors[5])


skip = 5
def update(num):
    ln1.set_data(QN[0:2, :num*skip])
    ln1.set_3d_properties(QN[2, :num*skip])

    ln2.set_data(QN[3*1:3*1+2, :num*skip])
    ln2.set_3d_properties(QN[3*1+2, :num*skip])

    ln3.set_data(QN[3*2:3*2+2, :num*skip])
    ln3.set_3d_properties(QN[3*2+2, :num*skip])

    ln4.set_data(QN[3*3:3*3+2, :num*skip])
    ln4.set_3d_properties(QN[3*3+2, :num*skip])

    ln5.set_data(QN[3*4:3*4+2, :num*skip])
    ln5.set_3d_properties(QN[3*4+2, :num*skip])

    ln6.set_data(QN[3*5:3*5+2, :num*skip])
    ln6.set_3d_properties(QN[3*5+2, :num*skip])


# Setting the axes properties


ani = animation.FuncAnimation(fig, update, N_anim, interval=1/100, blit=False)
writer = animation.PillowWriter(fps=25)
ani.save('C:/Users/Audre/Documents/Ponts/2A/S1/MOPSI/clean_code/MOPSI/Figures/animated/explicit.gif', writer=writer)

plt.show()
