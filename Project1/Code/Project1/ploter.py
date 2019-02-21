import numpy as np
import matplotlib.pyplot as plt

filePath = 'C:\\Users\\jonro\\Documents\\Comphys\\project1\\project1\\project1\\'

def E_analytical(alpha,N,D,omega):
    return D*N*alpha/2 + N*D*omega**2 /(8*alpha)

def analyticalVSnumerical(N,D):
    filName = filePath+ 'BfMcDim'+str(D) + 'Npart'+str(N) + 'Iter1000000.txt'
    list = np.loadtxt(filName)
    alpha = list[:,0]
    E_numerical = list[:,1]
    E_analytic = E_analytical(alpha,N,D,1)
    #variance = list[:,2]
    plt.figure('BruteForce'+'Dim' +str(D) +'N'+str(N))
    plt.plot(alpha,E_numerical,'ro',label='Numerical')
    #plt.plot(alpha,variance,'bo', label='Variance')
    plt.plot(alpha, E_analytic,label='Analytical')
    plt.legend()
    plt.show()
analyticalVSnumerical(10,1)
analyticalVSnumerical(10,2)
analyticalVSnumerical(10,3)