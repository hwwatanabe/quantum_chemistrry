from numpy import pi as PI
import math
import matplotlib.pyplot as plt
import numpy as np
"""
todo: 
structure opt
mixing beta
ab init md
"""

def main():

#    ### default ###
#    IOP = 2
#    N =  3
#    R = 1.4632
#    zeta1 = 2.0925
#    zeta2 = 1.24
#    ZA = 2 
#    ZB = 1
#
#    HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB)
#    ################
    IOP = 2
    N = 3 
    R = 1.4632
    zeta1 = 2.0925 
    zeta2 = 1.24
    ZA = 2 
    ZB = 1 

    HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB)

    plot_R_vs_E(IOP, N, zeta1, zeta2, ZA, ZB)


def plot_R_vs_E(IOP, N, zeta1, zeta2, ZA, ZB):

    R = 100
    energy_inf = HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB)

    energies = []
    Rs = np.linspace(0.001, 5, 500)
    for R in Rs:
        energy = HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB)
        energies.append(energy-energy_inf)
    plt.plot(Rs, energies)
    plt.ylim(-1,1)
    plt.axhline(0, ls="--", c="black")
    plt.show()


def HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB):

    ### calc integral ###

    R2 = R**2

    expon = np.zeros((3,3))
    coef = np.zeros((3,3))
    D1 = np.zeros(3)
    A1 = np.zeros(3)
    D2 = np.zeros(3)
    A2 = np.zeros(3)


    coef = [
            [1, 0, 0],
            [0.678914, 0.430129, 0],
            [0.444635, 0.535328, 0.154329]
            ]

    expon = [
            [0.270950, 0, 0],
            [0.151623, 0.851819, 0],
            [0.109818, 0.405771, 2.22766]
            ]

    for i in range(N):
        A1[i] = expon[N-1][i]*(zeta1**2)
        D1[i] = coef[N-1][i]*((2.0*A1[i]/PI)**0.75) 
        A2[i] = expon[N-1][i]*(zeta2**2)
        D2[i] = coef[N-1][i]*((2.0*A2[i]/PI)**0.75) 

    S12  = 0
    T11  = 0
    T12  = 0
    T22  = 0 
    V11A = 0
    V12A = 0
    V22A = 0
    V11B = 0
    V12B = 0
    V22B = 0
    V1111 = 0
    V2111 = 0 
    V2121 = 0 
    V2211 = 0
    V2221 = 0
    V2222 = 0

    for i in range(N):
        for j in range(N):
            RAP = A2[j]*R/(A1[i]+A2[j])
            RAP2 = RAP**2 
            RBP2 = (R-RAP)**2

            S12 += SX(A1[i], A2[j], R2)*D1[i]*D2[j]

            T11 += T(A1[i], A1[j], 0.0)*D1[i]*D1[j]
            T12 += T(A1[i], A2[j], R2 )*D1[i]*D2[j]
            T22 += T(A2[i], A2[j], 0.0)*D2[i]*D2[j]

            V11A += V(A1[i], A1[j], 0   , 0 , ZA)*D1[i]*D1[j]
            V12A += V(A1[i], A2[j], R2, RAP2, ZA)*D1[i]*D2[j]
            V22A += V(A2[i], A2[j], 0 , R2,   ZA)*D2[i]*D2[j]

            V11B += V(A1[i], A1[j], 0   , R2, ZB)*D1[i]*D1[j]
            V12B += V(A1[i], A2[j], R2, RBP2, ZB)*D1[i]*D2[j]
            V22B += V(A2[i], A2[j], 0 , 0 ,   ZB)*D2[i]*D2[j]

    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    RAP = A2[i]*R/(A2[i]+A1[j])
                    RBP = R - RAP
                    RAQ = A2[k]*R/(A2[k]+A1[l])
                    RBQ = R - RAQ
                    RPQ = RAP - RAQ

                    RAP2 = RAP**2
                    RBP2 = RBP**2
                    RAQ2 = RAQ**2
                    RBQ2 = RBQ**2
                    RPQ2 = RPQ**2

                    V1111 += twoe(A1[i], A1[j], A1[k], A1[l], 0, 0, 0    )*D1[i]*D1[j]*D1[k]*D1[l]
                    V2111 += twoe(A2[i], A1[j], A1[k], A1[l], R2, 0, RAP2)*D2[i]*D1[j]*D1[k]*D1[l]
                    V2121 += twoe(A2[i], A1[j], A2[k], A1[l], R2, R2,RPQ2)*D2[i]*D1[j]*D2[k]*D1[l]
                    V2211 += twoe(A2[i], A2[j], A1[k], A1[l], 0 , 0,R2   )*D2[i]*D2[j]*D1[k]*D1[l]
                    V2221 += twoe(A2[i], A2[j], A2[k], A1[l], 0 , R2,RBQ2)*D2[i]*D2[j]*D2[k]*D1[l]
                    V2222 += twoe(A2[i], A2[j], A2[k], A2[l], 0 , 0,0    )*D2[i]*D2[j]*D2[k]*D2[l]



    ### colect ###
    H = [
            [T11+V11A+V11B,T12+V12A+V12B],
            [T12+V12A+V12B,T22+V22A+V22B]
            ]

    H = np.array(H)


    S = [
            [1, S12],
            [S12, 1]
            ]

    S = np.array(S)

    X = [
            [1/math.sqrt(2*(1+S12)), 1/math.sqrt(2*(1-S12))],
            [1/math.sqrt(2*(1+S12)),-1/math.sqrt(2*(1-S12))]
            ]

    X = np.array(X)


    TT = np.zeros((2,2,2,2))

    TT[0][0][0][0] = V1111
    TT[1][0][0][0] = V2111
    TT[0][1][0][0] = V2111
    TT[0][0][1][0] = V2111
    TT[0][0][0][1] = V2111
    TT[1][0][1][0] = V2121
    TT[0][1][1][0] = V2121
    TT[1][0][0][1] = V2121
    TT[0][1][0][1] = V2121
    TT[1][1][0][0] = V2211
    TT[0][0][1][1] = V2211
    TT[1][1][1][0] = V2221
    TT[1][1][0][1] = V2221
    TT[1][0][1][1] = V2221
    TT[0][1][1][1] = V2221
    TT[1][1][1][1] = V2222
    
    if IOP != 0:
        print("aaa")



    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    print("{}{}{}{}: {}".format(i+1,j+1,k+1,l+1,TT[i,j,k,l]))


    ### scf ###
    crit = 1e-4
    maxit = 250


    niter = 0
    P = np.zeros((2,2))



# scf loop
    while ( niter <= maxit ):
        niter += 1

        G = np.zeros((2,2))
        for i in range(2):
            for j in range(2):
                G[i][j] = 0
                for k in range(2):
                    for l in range(2):
                        G[i][j] = G[i][j]+P[k][l]*(TT[i][j][k][l] - 0.5*TT[i][l][k][j])
    
    
        F = np.zeros((2,2))
        F = H + G
    
        energy = 0
        energy = 0.5*np.sum(P*(H+F))
    
    
        G = F@X

        FPRIME = X.T@G
        CPRIME, E = diag(FPRIME)
        C = X@CPRIME
#        print(C)
#        plot_phi(IOP, N, R, A1, D1, A2, D2, C)
    
        OLDP = P.copy()
        for i in range(2):
            for j in range(2):
                P[i,j] = 0
                for k in range(1):
                    P[i,j] += 2*C[i,k]*C[j,k]
    
        
        delta = 0
        delta += np.sum((P - OLDP)**2)
        delta = math.sqrt(delta/4)
        
    

        print("iter: {}, energy: {}, delta: {}".format(niter, energy, delta))
        OLDP = P@S
#        print(OLDP)

        if delta < crit:
            print("required accuracy has been reached, break scf loop")
            break
        
    if niter > maxit:
        print("not converged")
        return 0 
#        exit()

    energy_tot = energy + ZA*ZB/R

    print("###### final result #########")
    print("wave function:\n{}".format(C))
    print("eigen value:\n{}".format(E))
    print("Mulliken electron density:\n{}".format(P@S))

    C_temp, E_temp = diag(S)
    S_sqrt = C_temp@(E_temp**0.5)@C_temp.T
    PPRIME = S_sqrt@P@S_sqrt
    print("lowdin electron density:\n{}".format(PPRIME))
    print("electron energy: {}".format(energy))
    print("total energy: {}".format(energy_tot))



    return energy_tot 

def diag(F):
    if abs(F[0,0] - F[1,1]) < 1e-20:
        theta = PI/4
    else:
        theta = 0.5*math.atan(2*F[0,1]/(F[0,0]-F[1,1]))
#        theta = 0.5*math.atan2(2*F[0,1],(F[0,0]-F[1,1]))


    C = [
            [math.cos(theta), math.sin(theta)],
            [math.sin(theta), -math.cos(theta)]
        ]
    C = np.array(C)


    E = np.zeros((2,2))
    E[0,0] = F[0,0]*math.cos(theta)**2+F[1,1]*math.sin(theta)**2 + F[0,1]*math.sin(2*theta)
    E[1,1] = F[1,1]*math.cos(theta)**2+F[0,0]*math.sin(theta)**2 - F[0,1]*math.sin(2*theta)

    if E[1,1] < E[0,0]:
        temp = E[1,1]
        E[1,1] = E[0,0]
        E[0,0] = temp
        temp = C[0,1]
        C[0,1] = C[0,0]
        C[0,0] = temp
        temp = C[1,1]
        C[1,1] = C[1,0]
        C[1,0] = temp
    return [C, E]

def plot_phi(IOP, N, R, A1, D1, A2, D2, C):
    NN = 1000
    x = np.linspace(-R, 4*R, NN)
    phi1 = np.zeros(NN)
    phi2 = np.zeros(NN)
    for i in range(N):
        phi1 += D1[i]*np.exp(-1*A1[i]*x**2)
        phi2 += D2[i]*np.exp(-1*A2[i]*(x-R)**2)

    temp1 = phi1
    temp2 = phi2

    phi1 = C[0][0]*temp1 + C[1][0]*temp2
    phi2 = C[0][1]*temp1 + C[1][1]*temp2

    plt.plot(x, phi1, label="phi1_bonding")
    plt.plot(x, phi2, label="phi2_nonbonding")
#    plt.plot(x, phi1**2, "--", label="phi1_bonding")
#    plt.plot(x, phi2**2, "--", label="phi2_nonbonding")
    plt.legend()
    plt.ylim(-2,2)
    plt.scatter(R,0)
    plt.scatter(0,0)
    plt.axhline(0, ls="--", c="black")
    plt.show()

def F0(x):
    if not (x < 10**(-6)):
        temp = math.sqrt(PI/x)*math.erf(math.sqrt(x))/2
    else:
        temp = 1 - x/3

    return temp


def SX(A, B, RAB2):
    return (PI/(A+B))**1.5*math.exp(-1*A*B*RAB2/(A+B))


def T(A, B, RAB2):
    return A*B/(A+B)*(3-2*A*B*RAB2/(A+B))*(PI/(A+B))**1.5*math.exp(-1*A*B*RAB2/(A+B))


def V(A, B, RAB2, RCP2, ZC):
    temp = 2*PI/(A+B)*F0((A+B)*RCP2)*math.exp(-A*B*RAB2/(A+B))
    temp = -1*temp*ZC
    return temp


def twoe(A, B, C, D, RAB2, RCD2, RPQ2):
    temp = 2*(PI**2.5)/((A+B)*(C+D)*math.sqrt(A+B+C+D))*F0((A+B)*(C+D)*RPQ2/(A+B+C+D))*math.exp(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D))
    return temp


main()
