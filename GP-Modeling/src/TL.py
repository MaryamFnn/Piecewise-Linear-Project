# Lenght of transmission like a numpy array
#change the input signal to step signal
# what if I have pul parameters in matrix (four port network)


import numpy as np
import matplotlib.pyplot as plt
import sys

file = sys.argv[1]
#file = "example.txt"
Lzmin, Lzmax, num, L, C, R, G, Rg, Rl, Cl, A, th, tr, t0, nz, dt, nt, sen, snap = np.loadtxt(file,dtype='double')
Lens = np.linspace(Lzmin,Lzmax,num)
for Len in Lens:

    nz = int(nz)
    nt = int(nt)
    #L = Rc/v
    #C = 1/(Rc*v)
    th = th-tr
    V = np.zeros(2,nz+1)
    I = np.zeros(nz)
    dz = Lz/nz
    sen = np.argmin(np.fabs(np.arange(nz+1)*dz-sen))

    ciself = np.linalg.inv((L/dt+R/2)) * (L/dt-R/2)
    cvself = (C/dt-G/2)/(C/dt+G/2)
    ci = 1/((L/dt+R/2)*dz)
    cv = 1/((C/dt+G/2)*dz)

    cv0s = (C*dz/(2.0*dt)-1.0/(2.0*Rg))/(C*dz/(2.0*dt)+1.0/(2.0*Rg))
    cv0src = (1.0/Rg)/(C*dz/(2.0*dt)+1.0/(2.0*Rg))

    cv0i = 1.0/(C*dz/(2.0*dt)+1.0/(2.0*Rg))
    cvls = (1.0-(2.0*dt/(C*dz))*(1/(2.0*Rl) - Cl/dt))/(1.0+(2.0*dt/(C*dz))*(1/(2.0*Rl) + Cl/dt))
    cvli = (2.0*dt/(C*dz))/(1.0+(2.0*dt/(C*dz))*(1/(2.0*Rl) + Cl/dt))

    tt = np.arange(nt)*dt
    condition = [np.logical_and((tt-t0)<tr,tt-t0>0), np.logical_and(tt-t0>=tr,tt-t0<tr+th), np.logical_and(tt-t0>=tr+th,tt-t0<2*tr+th)]
    value = [lambda x: A*(x-t0)/tr, A, lambda x: A*(-(x-t0)/tr+2+th/tr), 0]

    vs = np.piecewise(tt, condition, value)
    vs = np.pad(0.5*(vs[1:]+vs[:-1]),(0,1))


    sens = np.zeros(nt)
    for it in range(nt):
        if (it==snap):
            plt.plot(np.arange(nz+1)*dz,V)
            #plt.plot((0.5+np.arange(nz))*dz,I*Rc)
            plt.ylim(-1,1)        
            plt.xlabel("z-position [m]")
            plt.ylabel("voltage [V]")
            plt.title("Voltage at t=%2.3f ns" %(1.0e9*it*dt))
            plt.savefig(file[:-4]+"time.png")
            #plt.show()
            plt.close()
        I = ciself*I-ci*(V[1:]-V[:-1])
        V[0] = cv0s*V[0] + cv0src*vs[it] - cv0i*I[0]
        V[1:-1] = cvself*V[1:-1] - cv*(I[1:]-I[:-1])
        V[-1] = cvls*V[-1] + cvli*I[-1]    
        sens[it]=1.0*V[0,-1]
        sens2[it]=1.0*V[1,-1]
        
    #plt.plot((np.arange(nt))*dt,vs)       
    plt.plot((np.arange(0,nt))*dt,sens)
    plt.xlabel("time [s]")
    plt.ylabel("Voltage [V]")
    plt.title("Voltage at z = %3.3em" %(np.arange(nz+1)[sen]*dz))
    plt.savefig(file[:-4]+"position.png")
    #plt.show()
    plt.close()
