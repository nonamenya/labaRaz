#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
e = -4.8e-10
ep = 4.8e-10
m = 1.67e-24
c = 3e10
pi = math.pi

n0 = 10e9
d = 1
B = 10000
omega = (e*B)/(m*c)
a = 8.86e15

t = 0
z = -0.5
v = 0
dt = 1e-11
E = m*v**2/2

u = v/c
tau = omega*t
ksi = z*omega/c
dtau = omega*dt

t_list = []
z_list = []
v_list = []
E_list = []

def E_field(z):
    return 4*pi*e*d*n0*z
print(e*E_field(0.5)/m)
print(math.sqrt(200/a))
while t < math.sqrt(200/a):
    z_list.append(z)
    v_list.append(v)
    t_list.append(t)
    E_list.append(E)

    tau = tau + dtau
    if z >= -d/2 and z <= d/2:
        u = u + (-E_field(z)/(B) - (a*m)/(e*B))*dtau/2
        ksi = ksi + u*dtau
        u = u + (-E_field(z)/(B) - (a*m)/(e*B))*dtau/2
    if z < -d/2:
        u = u + (E_field(-d/2)/(B) - (a*m)/(e*B))*dtau/2
        ksi = ksi + u*dtau
        u = u + (E_field(-d/2)/(B) - (a*m)/(e*B))*dtau/2
    if z >d/2:
        u = u + (E_field(d/2)/(B) - (a*m)/(e*B))*dtau/2
        ksi = ksi + u*dtau
        u = u + (E_field(d/2)/(B) - (a*m)/(e*B))*dtau/2
    v = c*u
    t = tau/omega
    z = ksi*c/omega
    E = m*v**2/(2*1.6e-12)  

with plt.style.context('bmh'):
    fig, axes = plt.subplots(4, 1, figsize = (10, 10))
    axes[0].plot(t_list, z_list)
    axes[0].set_xlabel('t, s')
    axes[0].set_ylabel('z, cm')

    

    axes[2].plot(z_list, E_list)
    axes[2].set_xlabel('z, cm')
    axes[2].set_ylabel('E, eV')
    
    axes[3].plot(t_list, E_list)
    axes[3].set_xlabel("t, s")
    axes[3].set_ylabel("E, eV")
    
    plt.savefig('result.png')
    plt.grid(True)
    plt.show()


# In[ ]:




