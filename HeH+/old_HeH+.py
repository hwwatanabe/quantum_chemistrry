import numpy as np
import math

def main():
    IOP = 2
    N = 3
    R = 1.4632
    zeta1 = 2.0925
    zeta2 = 1.24
    ZA = 2
    ZB = 1

    HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB)


def HFcalc(IOP, N, R, zeta1, zeta2, ZA, ZB):
    pass

    if IOP != 0:
        print("STO-{}G for atomic numbers {} and {}".format(N, ZA, ZB))

    intgrl(IOP, N, R, zeta1, zeta2, ZA, ZB)
    colect(IOP, N, R, zeta1, zeta2, ZA, ZB)
    scf(IOP, N, R, zeta1, zeta2, ZA, ZB)


def intgrl(IOP, N, R, zeta1, zeta2, ZA, ZB):
    R2 = R**2
    expon = np.zeros((3,3))
    coef = np.zeros((3,3))
    D1 = np.zeros(3)
    A1 = np.zeros(3)
    D2 = np.zeros(3)
    A2 = np.zeros(3)

    expon = [
            [1, 0, 0],
            [0.678914, 0.430129, 0],
            [0.444635, 0.535328, 0.154329]
            ]

    coef = [
            [0.270950, 0, 0],
            [0.151623, 0.851819, 0],
            [0.109818, 0.405771, 2.22766]
            ]


    for i in range(N):
        A1[i] = expon(i,N)*(zeta1**2)
        D1[i] = coef(i,N)*((2.0*A1(i)/np.pi)**0.75) 
        A2[i] = expon(i,N)*(zeta2**2)
        D2[i] = coef(i,N)*((2.0*A2(i)/np.pi)**0.75) 

    for i in range(N):
        for j in range(N):
            RAP = A2[j]*R/(A1[i]+A2[j])
            RAP2 = RAP**2 
            RBP2 = (R-RAP)**2

            S12 = S(A[i], A2[j], R2)*D1[i]*D2[j]

            T11 = T(A1[i], A1[j], 0.0)*D1[i]*D1[j]
            T12 = T(A1[i], A2[j], R2 )*D1[i]*D2[j]
            T22 = T(A2[i], A2[j], 0.0)*D2[i]*D2[j]

            V11A = V(A1[i], A1[j], 0   , 0 , ZA)*D1[i]*D1[j]
            V12A = V(A1[i], A2[j], R2, RAP2, ZA)*D1[i]*D2[j]
            V22A = V(A1[i], A2[j], 0 , R2,   ZA)*D2[i]*D2[j]

            V11B = V(A1[i], A1[j], 0   , 0 , ZB)*D1[i]*D1[j]
            V12B = V(A1[i], A2[j], R2, RBP2, ZB)*D1[i]*D2[j]
            V22B = V(A1[i], A2[j], 0 , R2,   ZB)*D2[i]*D2[j]
    
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    RAP = A2[i]*R/(A2[i]+A1[j])
                    RBP = R - RAP
                    RAQ = A2[k]*R/(A2[k]+A1[l])
                    RBQ = R - RAQ

                    RAP2 = RAP**2
                    RBP2 = RBP**2
                    RAQ2 = RAQ**2
                    RPQ2 = RPQ**2

                    V1111 = twoe(A1[i], A1[j], A1[k], A1[l], 0, 0, 0    )*D1[i]*D1[j]*D1[k]*D1[l]
                    V2111 = twoe(A2[i], A1[j], A1[k], A1[l], R2, 0, RAP2)*D2[i]*D1[j]*D1[k]*D1[l]
                    V2121 = twoe(A2[i], A1[j], A2[k], A1[l], R2, R2,RPQ2)*D2[i]*D1[j]*D2[k]*D1[l]
                    V2211 = twoe(A2[i], A2[j], A1[k], A1[l], 0 , 0,R2   )*D2[i]*D2[j]*D1[k]*D1[l]
                    V2221 = twoe(A2[i], A2[j], A2[k], A1[l], 0 , R2,RBQ2)*D2[i]*D2[j]*D1[k]*D1[l]
                    V2222 = twoe(A2[i], A2[j], A2[k], A2[l], 0 , 0,0    )*D2[i]*D2[j]*D2[k]*D2[l]


    
    if IOP != 0:
        print("aaaaaaa")


def F0(x):
    if not (x < 10**(-6)):
        temp = math.sqrt(np.pi/x)*math.erf(math.sqrt(x))/2
    else:
        temp = 1 - x/3

    return temp:


def S(A, B, RAB2):
    return np.pi/(A+B)**1.5*math.exp(-1*A*B*RAB2/(A+B))


def T(A, B, RAB2):
    return A*B/(A+B)*(3e-2*A*B*RAB2/(A+B))*(np.pi/(A+B))**1.5*math.exp(-1*A*B*RAB2/(A+B))


def V(A, B, RAB2, RCP2, ZC):
    temp = 2*np.pi/(A+B)*F0((A+B)*RCP2)*math.exp(-A*B*RAB2/(A+B))
    temp = -1*temp*ZC
    return temp

def twoe(A, B, C, D, RAB2, RCD2, RPQ2):
    temp = 2*(np.pi**2.5)/((A+B)*(C+D)*math.sqrt(A+B+C+D))*F0((A+B)*(C+D)*RPQ2/(A+B+C+D))*math.exp(-A*B*RAB2/(A+B)-C*D*RC2/(C+D))
    return temp


def colect(IOP, N, R, zeta1, zeta2, ZA, ZB):
    H = [
            [T11+V11A+V11B,T12+V12A+V12B],
            [T12+V12A+V12B,T22+V22A+V22B]
            ]

    S = [
            [1, S12],
            [S12, 1]
            ]


    X = [
            [1/math.sqrt(2*(1+S12)), 1/math.sqrt(2*(1+S12))],
            [1/math.sqrt(2*(1-S12)), -1/math.sqrt(2*(1-S12))]
            ]

    XT = [
            [X[0][0], X[1][0]],
            [X[0][1], X[1][1]]
            ]

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
    TT[0][1][0][1] = V2211
    TT[0][0][1][1] = V2211
    TT[1][1][1][0] = V2221
    TT[1][1][0][1] = V2221
    TT[1][0][1][1] = V2221
    TT[0][1][1][1] = V2221
    TT[1][1][1][1] = V2222
    
    if IOP != 0:
        matout(S, 2,2,2,2, 4HS)
        matout(X, 2,2,2,2, 4HX)
        matout(H, 2,2,2,2, 4HH)



    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    print("aaa")

        
def scf(IOP, N, R, zeta1, zeta2, ZA, ZB):
    crit = 1e-4
    maxit = 25

    niter = 0

    p = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            P[i][j] = 0

    if IOP == 2:
        matout(P, 2,2,2,2, 4HP)

    niter += 1

    if IOP == 2:
        print("aaa")

    formg()

    if IOP == 2:
        matout(G, 2,2,2,2, 4HG)

    for i in range(2):
        for j in range(2):
            F[i][j] = H[i][j] + G[i][j]

    en = 0
    for i in range(2):
        for j in range(2):
        en += 0.5*P[i][j]*(H[i][j] + F[i][j])

    if IOP ==2 :
        matout(F,2,2,2,2, 4HF)
        print("aaa")

    mult(F,X, G, 2,2)
    mult(XT, G, FPRIME, 2, 2)
    diag(FPRIME, CPRIME, E)
    mult(X, CPRIME, C, 2,2)

    OLDP = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            OLDP[i][j] = P[i][j]
            P[i][j] = 0

            for k in range(1):
                P[i][j] += 2*C[i][k]*C[J][k]

    if IOP == 2:
        matout(FPRIME, 2,2,2,2, 4HF)
        matout(CPRIME, 2,2,2,2, 4HF)


            







