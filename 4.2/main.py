import math
class Prog:
    def __init__(self, N1, N2):
        self.N1 = N1
        self.N2 = N2
        self.res = [.0] * (N2 + 1)
        self.h = 1/N1
        self.tao = 1/N2
        for i in range(N2+1):
            self.res[i] = [.0] * (N1 + 1)

    def sigma(self):
        return ((1/4)-(self.h*self.h)/(12*self.tao*self.tao))
    def f(self, x, t):
        return ((x*x-t*t)*math.exp(-1*x*t))
    def f__(self, x,t):
        return ((2-4*t*x+t*t*x*x-t*t*t*t) * math.exp(-1 * x * t))
    def u0(self, x):
        return 1
    def u1(self, x):
        return (-x)
    def miu0(self, t):
        return(1)
    def miu1(self, t):
        return math.exp(-t)
    def phi(self, x, t):
        return (self.f(x,t)+((self.h*self.h)/12)*self.f__(x,t))
    def U_1(self, x):
        return (self.u1(x)+(self.tao/2)*self.f(x,0))

    def A(self):
        return (self.sigma()/(self.h*self.h))
    def B(self):
        return (self.sigma() / (self.h * self.h))
    def C(self):
        return ((2*self.sigma() / (self.h * self.h))+(1/(self.tao*self.tao)))
    def F(self, i, j):
        F = self.phi(self.h*i,self.tao*j)+(1/(self.tao*self.tao))*(2*self.res[j][i]-self.res[j-1][i])+(1-2*self.sigma())*(1/(self.h*self.h))*(self.res[j][i+1]-2*self.res[j][i]+self.res[j][i-1])+self.sigma()*(1/(self.h*self.h))*(self.res[j-1][i+1]-2*self.res[j-1][i]+self.res[j-1][i-1])
        return F
#####Реализация метода прогонки
    def Progonka(self):

        for j in range(self.N2+1):
            t_j = self.tao*j
            self.res[j][0] = self.miu0(t_j)
            self.res[j][self.N1] = self.miu1(t_j)
        for i in range(1, self.N1):
            x_i = self.h*i
            self.res[0][i] = self.u0(x_i)
            self.res[1][i] = self.u0(x_i)+self.U_1(x_i)*self.tao

        alpha = [.0]*(self.N1+1)
        beta = [.0]*(self.N1+1)
        for j in range(1, self.N2):
            alpha[1] = 0
            beta[1] = self.miu0(j*self.tao)
            for i in range(1, self.N1):
                alpha[i+1] = self.B()/(self.C()-self.A()*alpha[i])
                beta[i+1] = (self.F(i,j)+beta[i]*self.A())/(self.C()-self.A()*alpha[i])

            for i in range((self.N1-1), 0, -1):
                self.res[j+1][i] = alpha[i+1]*self.res[j+1][i+1]+beta[i+1]
        return (self.res)

def print_res(res, N1, N2):
    for j in range(N2+1):
        for i in range(N1+1):
            print(f"{(res[j][i]):.4f}", end = " ")
        print("")
def print_pogr(res, res_t, N1, N2):
    for j in range(N2+1):
        for i in range(N1+1):
            print(f"{abs(res[j][i]-res_t[j][i]):.1E}", end = " ")
        print("")
def Test20():
    p = Prog(20, 20)
    res = p.Progonka()
    print("Решение при количестве разбиений 20:")
    print_res(res, 20, 20)
    p_t = Prog(400, 400)
    res_t1000 = p_t.Progonka()
    res_t = [.0] * 401
    for i in range(400 + 1):
        res_t[i] = [.0] * (400 + 1)
    for j in range(21):
        for i in range(21):
            res_t[j][i] = res_t1000[j * 20][i * 20]
    print("Решение при количестве разбиений 400:")
    print_res(res_t, 20, 20)
    print("Матрица невязок полученных решений:")
    print_pogr(res, res_t, 20, 20)
def Test10():
    p = Prog(10, 10)
    res = p.Progonka()
    print("Решение при количестве разбиений 10:")
    print_res(res, 10, 10)
    p_t = Prog(100, 100)
    res_t1000 = p_t.Progonka()
    res_t = [.0] * 101
    for i in range(100 + 1):
        res_t[i] = [.0] * (100 + 1)
    for j in range(11):
        for i in range(11):
            res_t[j][i] = res_t1000[j * 10][i * 10]
    print("Решение при количестве разбиений 100:")
    print_res(res_t, 10, 10)
    print("Матрица невязок полученных решений:")
    print_pogr(res, res_t, 10, 10)
###################main
Test20()
Test10()