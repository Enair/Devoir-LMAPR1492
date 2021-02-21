import matplotlib.pyplot as plt
import math as m
import numpy as np

b = 3*int(m.sqrt(58))
x1 = np.arange(-2*np.pi, -3*np.pi/2, 0.01)
x2 = np.arange(-3*np.pi/2, -np.pi/2, 0.01)
x3 = np.arange(-np.pi/2, np.pi/2, 0.01)
x4 = np.arange(np.pi/2, 3*np.pi/2, 0.01)
x5 = np.arange(3*np.pi/2, 2*np.pi, 0.01)



k1 = 0.5
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
ax.plot(x1,np.real(s1), color = 'blue', label = "real")
ax.plot(x2,np.real(s2), color = 'blue', label = "real")
ax.plot(x3,np.real(s3), color = 'blue', label = "real")
ax.plot(x4,np.real(s4), color = 'blue', label = "real")
ax.plot(x5,np.real(s5), color = 'blue', label = "real")

ax.plot(x1,np.imag(s1), color = 'green')
ax.plot(x2,np.imag(s2), color = 'green')
ax.plot(x3,np.imag(s3), color = 'green')
ax.plot(x4,np.imag(s4), color = 'green')
ax.plot(x5,np.imag(s5), color = 'green')

ax.set(xlabel = 'x', ylabel = 'phi_k(r)', title= 'fonction de bloch pour k = (1/2,0), pour r allant de [-2,0] à [2,0]')
ax.grid()
fig.savefig("phi 4")
plt.show()

