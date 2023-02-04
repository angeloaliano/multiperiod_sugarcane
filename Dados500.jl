using Printf, Random,JuMP, Gurobi, Printf
using FileIO, JLD2, LinearAlgebra
using Statistics
rng = MersenneTwister(2022);

ka = 400
ke = 100

Va = 1:20
Ve = 21:25
V = union(Va,Ve)
Ja = 1:400
Je = 401:500
J = union(Ja,Je)

H = [1,2,3,4,9,10]
M = [4,5,6,7,8,9,10,11]

PSC = [132.8, 121.93, 137.2, 89.29, 121.7, 140.6, 129.10, 130.90, 100, 112.30, 165.00, 117.80, 148.20, 113.00, 123.10, 142.40, 132.12, 112.80, 126.70, 136,0,0,0,0,0]
PEC = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,236, 215, 200, 181, 179]

SSP = [0.145,0.14,0.15,0.148,0.148,0.15,0.143,0.1354,0.1584,0.15,0.135,0.145,0.14,0.156,0.1332,0.157,0.14,0.149,0.132,0.15,0,0,0,0,0]
SEP = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.08,0.06,0.07,0.09,0.08]

FSP = [0.124,0.1238,0.115,0.123,0.113,0.129,0.122,0.12,0.1234,0.1238,0.115,0.124,0.12,0.124,0.12,0.1293,0.1238,0.129,0.1274,0.1116,0,0,0,0,0]
FEP = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.23,0.24,0.32,0.22,0.25]

L = [20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
     25, 30,20, 25, 30, 20, 25,
     20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
     25, 30,20, 25, 30, 20, 25,20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
     25, 30,20, 25, 30, 20, 25,
     20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
     25, 30,20, 25, 30, 20, 25,20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
     30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
      25, 30,20, 25, 30, 20, 25,
      20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
      30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
      25, 30,20, 25, 30, 20, 25,20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
      20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
      30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
      25, 30,20, 25, 30, 20, 25,
      20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25,
      30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20,
      25, 30,20, 25, 30, 20, 25,20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25, 30, 20, 25]

L=L[1:500]
DS = [zeros(3,4);20000*ones(8,4);zeros(1,4)]
DF = [zeros(3,4);20000*ones(8,4);zeros(1,4)]
Cap = [zeros(3,4);270000*ones(8,4);zeros(1,4)]

D=[-3,-2,-1,0,1,2,3]
C = [1,2,3,4]

αa = 0.06
αe = 0.06

SPSH(i,d,c) = (-0.0243*d^2 + 1)*((1-αa)^(c-1))*SSP[i]*PSC[i]
SPEH(i,d,c) = (0.0041*d + 1)*((1+αe)^(c-1))*SEP[i]*PEC[i]
FPSH(i,d,c) = (-0.0243*d^2 + 1)*((1-αa)^(c-1))*FSP[i]*PSC[i]
FPEH(i,d,c) = (0.0041*d + 1)*((1+αe)^(c-1))*FEP[i]*PEC[i]

#SPSH(i,d,c) = (-0.0243*d^2 + 1)*((1-αa)^(c-1))*SSP[i]*PSC[i]
#SPEH(i,d,c) = (-0.0243*d^2 + 1)*((1-αa)^(c-1))*SEP[i]*PSC[i]
#FPSH(i,d,c) = (0.0041*d + 1)*((1+αe)^(c-1))*FSP[i]*PEC[i]
#FPEH(i,d,c) = (0.0041*d + 1)*((1+αe)^(c-1))*FEP[i]*PEC[i]


#Custo = [57.57, 60.30, 62.54, 65.92, 66.35,79, 83.37, 77.2, 102.9, 83.47, 76.28, 80.41, 79.71, 95.07, 87.94, 69.44, 85.24, 73.90, 87.58, 82.86, 86.81, 45.70, 87.69, 81.36, 77.84]

CA = zeros(ka+ke,4)
for j in 1:ka
    for c=1
        CA[j,c] = 94
    end
    for c in [2,3,4]
        CA[j,c] = 65
    end
end
CE = zeros(ka+ke,4)
for j in ka+1:ka+ke
    for c=1
        CE[j,c] = 94
    end
    for c in [2,3,4]
        CE[j,c] = 65
    end
end
