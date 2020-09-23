#!/usr/bin/env python
# coding: utf-8

# In[41]:


import math 
import matplotlib.pyplot as plt 
e = 4.8e-10 
q=-1*e
m = 9.1e-28 
c = 3e+10 
E = 0.33 
f = 2.45e+9 
En = 2
gamma = (En*1.6e-12/(m*c**2))+1 
fi = math.pi/6 
B0 = (2*math.pi*f*m*c)/e 
w = (q*B0)/(m*c)
t = 0 
tau = t*w
dtau = 2*math.pi/300 
T = 2*math.pi/w
g0 = (E*e)/(m*c*2*math.pi*f)
g00=g0
alphamax = 1.9*g0**(4/3)
alpha = 0.1*alphamax 
omega=0.1
beta=0.1
B = B0*(1+alpha*tau)
b = B/B0 
t_list = [] 
fi_list = [] 
gamma_list = [] 
def fun_gamma (gamma, fi, b): 
  return -g0*(1 - 1/(gamma**2))**0.5*math.cos(fi) + 0.5*alpha*(1 - 1/(gamma**2))*gamma/b 
def fun_fi(fi, gamma, b): 
  return (b-gamma)/gamma + g0*((gamma**2 - 1)**(-0.5))*math.sin(fi) 

def ruku4_gamma (dy, y, fi, b): 
  ku1 = fun_gamma (y, fi, b) 
  ku2 = fun_gamma (y + dtau*0.5*ku1, fi, b) 
  ku3 = fun_gamma (y + dtau*0.5*ku2, fi, b) 
  ku4 = fun_gamma (y + dtau*ku3, fi, b) 
  dy = dy + (ku1 + 2*ku2 + 2*ku3 + ku4)*dtau/6 
  return dy 
def ruku4_fi (dy, y, gamma, b): 
  ku1 = fun_fi (y, gamma, b) 
  ku2 = fun_fi (y + dtau*0.5*ku1, gamma, b) 
  ku3 = fun_fi (y + dtau*0.5*ku2, gamma, b) 
  ku4 = fun_fi (y + dtau*ku3, gamma, b) 
  dy = dy + (ku1 + 2*ku2 + 2*ku3 + ku4)*dtau/6 
  return dy 
while tau<5000*math.pi: 
 
  
  fi_list.append(fi) 
  gamma_list.append(En) 
  t_list.append (t) 
  tau = tau + dtau 
  B = B0*(1+alpha*tau)
  b = B/B0 
  gamma = ruku4_gamma (gamma, gamma, fi, b) 
  fi = ruku4_fi(fi, fi, gamma, b) 
  En = (m*c**2*(gamma-1))/(1.6e-12) 
  t = tau/w/T

plt.figure (figsize = (11,7))
plt.plot(t_list, fi_list) 
plt.savefig('result.png') 
plt.xlabel('tau (s)') 
plt.ylabel('En (Ev)') 
plt.grid(True) 
plt.show()


# In[38]:


import math 
import matplotlib.pyplot as plt 
e = 4.8e-10 
m = 9.1e-28 
c = 3e+10 
E = 0.33 
f = 2.45e+9 
En = 2000 
gamma = (En*1.6e-12/(m*c**2))+1 
fi = math.pi/6 
B0 = (2*math.pi*f*m*c)/e 
w = (q*B0)/(m*c)
t = 0 
tau = t*w
dtau = 2*math.pi
T = 2*math.pi/w
g0 = (E*e)/(m*c*2*math.pi*f)
alphamax = 1.9*g0**(4/3)
alpha = 0.1*alphamax 
omega=0.1
beta=0.1
B = B0*(1+alpha*tau+beta*math.sin(omega*tau)) 
b = B/B0 
t_list = [] 
fi_list = [] 
gamma_list = [] 
def fun_gamma (gamma, fi, b): 
  return -g0*(1 - 1/(gamma**2))**0.5*math.cos(fi) + 0.5*alpha*(1 - 1/(gamma**2))*gamma/b 
def fun_fi(fi, gamma, b): 
  return (b-gamma)/gamma + g0*((gamma**2 - 1)**(-0.5))*math.sin(fi) 

def ruku4_gamma (dy, y, fi, b): 
  ku1 = fun_gamma (y, fi, b) 
  ku2 = fun_gamma (y + dtau*0.5*ku1, fi, b) 
  ku3 = fun_gamma (y + dtau*0.5*ku2, fi, b) 
  ku4 = fun_gamma (y + dtau*ku3, fi, b) 
  dy = dy + (ku1 + 2*ku2 + 2*ku3 + ku4)*dtau/6 
  return dy 
def ruku4_fi (dy, y, gamma, b): 
  ku1 = fun_fi (y, gamma, b) 
  ku2 = fun_fi (y + dtau*0.5*ku1, gamma, b) 
  ku3 = fun_fi (y + dtau*0.5*ku2, gamma, b) 
  ku4 = fun_fi (y + dtau*ku3, gamma, b) 
  dy = dy + (ku1 + 2*ku2 + 2*ku3 + ku4)*dtau/6 
  return dy 
while  tau<25000*math.pi:
    
    fi_list.append(fi) 
    gamma_list.append(En) 
    t_list.append (t) 
    tau = tau + dtau 
    B = B0*(1+alpha*tau) 
    b = B/B0 
    gamma = ruku4_gamma (gamma, gamma, fi, b)
    fi = ruku4_fi(fi, fi, gamma, b) 
    En = (m*c**2*(gamma-1))/(1.6e-12) 
    t = tau/w/T/200
    
    

plt.figure (figsize = (11,7))
plt.plot(t_list, gamma_list) 
plt.savefig('result.png') 
plt.xlabel('tau (s)') 
plt.ylabel('En (Ev)') 
plt.grid(True) 
plt.show()


# In[42]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math 
import matplotlib.pyplot as plt 
e = 4.8e-10 
q = -1*e
m = 9.1e-28 
c = 3e+10 
E = 0.33 
f = 2.45e+9 
En = 2000 
gamma = (En*1.6e-12/(m*c**2))+1 
fi = math.pi/6 
B0 = (2*math.pi*f*m*c)/e 
w = (q*B0)/(m*c)
t = 0 
tau = t*w
dtau = 2*math.pi
T = 2*math.pi/w
g0 = (E*e)/(m*c*2*math.pi*f)
alphamax = 1.9*g0**(4/3)
alpha = 0.1*alphamax 
omega=0.1
beta=0.1
B = B0*(1+alpha*tau+beta*math.sin(omega*tau)) 
b = B/B0 
t_list = [] 
fi_list = [] 
gamma_list = [] 
def fun_gamma (gamma, fi, b): 
  return -g0*(1 - 1/(gamma**2))**0.5*math.cos(fi) + 0.5*alpha*(1 - 1/(gamma**2))*gamma/b 
def fun_fi(fi, gamma, b): 
  return (b-gamma)/gamma + g0*((gamma**2 - 1)**(-0.5))*math.sin(fi) 

def ruku4_gamma (dy, y, fi, b): 
  ku1 = fun_gamma (y, fi, b) 
  ku2 = fun_gamma (y + dtau*0.5*ku1, fi, b) 
  ku3 = fun_gamma (y + dtau*0.5*ku2, fi, b) 
  ku4 = fun_gamma (y + dtau*ku3, fi, b) 
  dy = dy + (ku1 + 2*ku2 + 2*ku3 + ku4)*dtau/6 
  return dy 
def ruku4_fi (dy, y, gamma, b): 
  ku1 = fun_fi (y, gamma, b) 
  ku2 = fun_fi (y + dtau*0.5*ku1, gamma, b) 
  ku3 = fun_fi (y + dtau*0.5*ku2, gamma, b) 
  ku4 = fun_fi (y + dtau*ku3, gamma, b) 
  dy = dy + (ku1 + 2*ku2 + 2*ku3 + ku4)*dtau/6 
  return dy 
while  B<4*B0:
    
    fi_list.append(fi) 
    gamma_list.append(En) 
    t_list.append (t) 
    tau = tau + dtau 
    B = B0*(1+alpha*tau) 
    b = B/B0
    if B>=2*B0:
        b = 1
 
    gamma = ruku4_gamma (gamma, gamma, fi, b)
    fi = ruku4_fi(fi, fi, gamma, b) 
    En = (m*c**2*(gamma-1))/(1.6e-12) 
    t = tau/w/T/200
    
    

plt.figure (figsize = (11,7))
plt.plot(t_list, fi_list) 
plt.savefig('result.png') 
plt.xlabel('tau (s)') 
plt.ylabel('En (Ev)') 
plt.grid(True) 
plt.show()


# In[ ]:


if B>2*B0 and g0>0.0:
        g0=0
    elif g0>1.0e-12:
        g0=g00*(1-(tau/250)**2)
    if B>4*B0:
        alpha = 0
        B0 = B

