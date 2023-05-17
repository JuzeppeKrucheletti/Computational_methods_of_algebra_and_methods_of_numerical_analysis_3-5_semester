import math

class Jacobi:
    def __init__(self, a, b, c, d, h1, h2, h1_m, h2_m):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.h1 = h1
        self.h2 = h2
        self.N1 = int((self.b - self.a)/self.h1)
        self.N2 = int((self.d - self.c)/self.h2)
        self.E = max(math.pow(self.h1, 3), math.pow(self.h2, 3))
        self.res = [.0]*(self.N2+1)
        for i in range(self.N2+1):
            self.res[i] = [.0]*(self.N1+1)
        self.F = [.0]*(self.N2+1)
        for i in range(self.N2+1):
            self.F[i] = [.0]*(self.N1+1)

        self.h1_m = h1_m
        self.h2_m = h2_m
        self.N1_m = int((self.b - self.a) / self.h1_m)
        self.N2_m = int((self.d - self.c) / self.h2_m)
        self.res_m = [.0] * (self.N2_m + 1)
        for i in range(self.N2_m + 1):
            self.res_m[i] = [.0] * (self.N1_m + 1)
        self.F_m = [.0] * (self.N2_m + 1)
        for i in range(self.N2_m + 1):
            self.F_m[i] = [.0] * (self.N1_m + 1)
        self.E_m = max(math.pow(self.h1_m, 3), math.pow(self.h2_m, 3))


    def psi1(self, y):
        return (math.sin(math.pi*y)*math.sin(math.pi*y))
    def psi2(self, y):
        return 0
    def psi3(self, x):
        return (math.cosh(x*x-3*x)-1)
    def psi4(self, x):
        return 0
    def f(self,x,y):
        return (math.cosh(x-y))
    def Solve_Jacobi(self):
        res1 = [.0] * (self.N2 + 1)
        for i in range(self.N2 + 1):
            res1[i] = [.0] * (self.N1 + 1)

        for j in range(self.N2+1):
            y_j = j*self.h2
            self.res[j][0] = self.psi1(y_j)
            res1[j][0] = self.psi1(y_j)
            self.res[j][self.N1] = self.psi2(y_j)
            res1[j][self.N1] = self.psi2(y_j)
        for i in range(self.N1+1):
            x_i = i*self.h1
            self.res[0][i] = self.psi3(x_i)
            res1[0][i] = self.psi3(x_i)
            self.res[self.N2][i] = self.psi4(x_i)
            res1[self.N2][i] = self.psi4(x_i)
        for j in range(self.N2+1):
            for i in range(self.N1+1):
                x_i = self.h1*i
                y_j = self.h2*j
                self.F[j][i] = self.f(x_i, y_j)

        self.count = 0

        while True:
            end = True
            self.count+=1

            for j in range(1,self.N2):
                for i in range(1, self.N1):
                    y_ij = (1/(2/(self.h1*self.h1)+2/(self.h2*self.h2)))*((res1[j][i+1]+res1[j][i-1])/(self.h1*self.h1)+(res1[j+1][i]+res1[j-1][i])/(self.h2*self.h2)+self.F[j][i])
                    self.res[j][i] = y_ij
            end = True
            for j in range(self.N2 + 1):
                for i in range(self.N1 + 1):
                    if(abs(self.res[j][i]-res1[j][i])>self.E):
                        end = False
                        break
            if(end):
                break
            else:
                for j in range(self.N2 + 1):
                    for i in range(self.N1 + 1):
                        res1[j][i] = self.res[j][i]
        return self.res
    def Solve_Jacobi_m(self):
        res1 = [.0] * (self.N2_m + 1)
        for i in range(self.N2_m + 1):
            res1[i] = [.0] * (self.N1_m + 1)

        for j in range(self.N2_m+1):
            y_j = j*self.h2_m
            self.res_m[j][0] = self.psi1(y_j)
            res1[j][0] = self.psi1(y_j)
            self.res_m[j][self.N1_m] = self.psi2(y_j)
            res1[j][self.N1_m] = self.psi2(y_j)
        for i in range(self.N1_m+1):
            x_i = i*self.h1_m
            self.res_m[0][i] = self.psi3(x_i)
            res1[0][i] = self.psi3(x_i)
            self.res_m[self.N2_m][i] = self.psi4(x_i)
            res1[self.N2_m][i] = self.psi4(x_i)
        for j in range(self.N2_m+1):
            for i in range(self.N1_m+1):
                x_i = self.h1_m*i
                y_j = self.h2_m*j
                self.F_m[j][i] = self.f(x_i, y_j)

        self.count_m = 0
        while True:
            end = True
            self.count_m += 1
            print(str(self.count_m))
            for j in range(1,self.N2_m):
                for i in range(1, self.N1_m):
                    y_ij = (1/(2/(self.h1_m*self.h1_m)+2/(self.h2_m*self.h2_m)))*((res1[j][i+1]+res1[j][i-1])/(self.h1_m*self.h1_m)+(res1[j+1][i]+res1[j-1][i])/(self.h2_m*self.h2_m)+self.F_m[j][i])
                    self.res_m[j][i] = y_ij
            end = True
            for j in range(self.N2_m + 1):
                for i in range(self.N1_m + 1):
                    if(abs(self.res_m[j][i]-res1[j][i])>self.E_m):
                        end = False
                        break
            if(end):
                break
            else:
                for j in range(self.N2_m + 1):
                    for i in range(self.N1_m + 1):
                        res1[j][i] = self.res_m[j][i]
        return self.res_m
    def print_res(self):
        print("Решение задачи дирихле методом Якоби")
        for j in range(self.N2 + 1):
            for i in range(self.N1 + 1):
                print(f"{self.res[j][i]:.3f}", end = " ")
            print("")
        print("Число итераций: "+str(self.count))
    def print_pogr(self):
        for j in range(self.N2 + 1):
            for i in range(self.N1 + 1):
                print(f"{self.res_m[j*int(self.N2_m/self.N2)][i*int(self.N1_m/self.N1)]:.3f}", end = " ")
            print("")
        for j in range(self.N2 + 1):
            for i in range(self.N1 + 1):
                print(f"{abs(self.res_m[j * int(self.N2_m / self.N2)][i * int(self.N1_m / self.N1)] - self.res[j][i]):.2E}", end=" ")
            print("")
        print("Число итераций для шагов h1_m и h2_m: " + str(self.count))


##########main
J = Jacobi(0,3,0,1,0.05,0.1, 0.005, 0.01)
J.Solve_Jacobi()
J.Solve_Jacobi_m()
J.print_res()
J.print_pogr()