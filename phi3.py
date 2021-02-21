import matplotlib.pyplot as plt
import math as m
import numpy as np

b = 3*int(m.sqrt(58))
x1 = np.arange(-2*np.pi, -3*np.pi/2, 0.01)
x2 = np.arange(-3*np.pi/2, -np.pi/2, 0.01)
x3 = np.arange(-np.pi/2, np.pi/2, 0.01)
x4 = np.arange(np.pi/2, 3*np.pi/2, 0.01)
x5 = np.arange(3*np.pi/2, 2*np.pi, 0.01)



k1 = 0.0
k2 = 0.0
r1 = [x1, x2, x3, x4, x5]
r2 = [0.0, 0.0, 0.0, 0.0, 0.0]

rnorm = [np.zeros(len(x1)),np.zeros(len(x2)),np.zeros(len(x3)),np.zeros(len(x4)),np.zeros(len(x5))]
for i in range(len(x1)):
    rnorm[0][i] = m.sqrt((x1[i]+np.pi*2)**2 + r2[0]**2)
for i in range(len(x2)):
    rnorm[1][i] = m.sqrt((x2[i]+np.pi)**2 + r2[0]**2)
for i in range(len(x3)):
    rnorm[2][i] = m.sqrt(x3[i]**2 + r2[0]**2)
for i in range(len(x4)):
    rnorm[3][i] = m.sqrt((x4[i]-np.pi)**2 + r2[0]**2)
for i in range(len(x5)):
    rnorm[4][i] = m.sqrt((x5[i]-np.pi*2)**2 + r2[0]**2)
    
    


s1 = np.exp(1j*(k1*r1[0] + k2*r2[0]))*np.exp(-b*(rnorm[0])**2)
s2 = np.exp(1j*(k1*r1[1] + k2*r2[1]))*np.exp(-b*(rnorm[1])**2)
s3 = np.exp(1j*(k1*r1[2] + k2*r2[2]))*np.exp(-b*(rnorm[2])**2)
s4 = np.exp(1j*(k1*r1[3] + k2*r2[3]))*np.exp(-b*(rnorm[3])**2)
s5 = np.exp(1j*(k1*r1[4] + k2*r2[4]))*np.exp(-b*(rnorm[4])**2)

fig, ax = plt.subplots()
ax.plot(np.real(s1),x1, color = 'blue', label = "real")
ax.plot(np.real(s2),x2, color = 'blue', label = "real")
ax.plot(np.real(s3),x3, color = 'blue', label = "real")
ax.plot(np.real(s4),x4, color = 'blue', label = "real")
ax.plot(np.real(s5),x5, color = 'blue', label = "real")



ax.set(xlabel = 'x', ylabel = 'phi_k(r)', title= 'fonction de bloch pour k = (0,0), pour r allant de [0,-2] à [0,2]')
ax.grid()
fig.savefig("phi 3")
plt.show()

