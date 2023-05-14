import math
################1
def K(x):
    return (1-x*x)
def K_(x):
    return (-2*x)
def F(x):
    return(math.cos(x)-2*x*math.sin(x))
def q(x):
    return(x*x)
def kappa0():
    return(1/2)
def g0():
    return(1/2)
def kappa1():
    return(1)
def g1():
    return (math.cos(1))
#######
def A(x, h):
    return ((-K_(x)/(2*h))+K(x)/(h*h))
def B(x, h):
    return ((K_(x) / (2 * h)) + K(x) / (h*h))
def C(x, h):
    return((2*K(x))/(h*h)+q(x))
def K1(h):
    return(K(0)/(h*(kappa0()-(h/2)*kappa0()*(K_(0)/K(0))+(h/2)*q(0))+K(0)))
def V1(h):
    return((h*(g0()-(h/2)*g0()*(K_(0)/K(0))+(h/2)*F(0)))/(h*(kappa0()-(h/2)*kappa0()*(K_(0)/K(0))+(h/2)*q(0))+K(0)))
def Progonka_1(N):
    h = 1/N
    Y = [.0]*(N+1)
    alpha = [.0]*(N+1)
    beta = [.0]*(N+1)
    alpha[1] = K1(h)
    beta[1] = V1(h)
    for i in range(1, N):
        x_i = i*h
        alpha[i+1] = B(x_i,h)/(C(x_i, h)-A(x_i, h)*alpha[i])
        beta[i+1] = (F(x_i)+beta[i]*A(x_i, h))/(C(x_i, h)-A(x_i, h)*alpha[i])
    Y[N] = g1()/kappa1()

    for i in range((N-1), -1, -1):
        Y[i] = alpha[i+1]*Y[i+1]+beta[i+1]
    return (Y)

def print_Y(Y):
    for i in range(len(Y)):
        print(str(Y[i]))
def print_pogr(Y, Y_t):
    for i in range(len(Y)):
        print(f"{abs(Y[i]-Y_t[i]):E}")
####################2
def a(x, h):
    return (K(x-(h/2)))
def phi(x,h):
    return (F(x))
def d(x,h):
    return (q(x))
def d0(h):
    return (q(h/4))
def phi0(h):
    return (F(h/4))

def A_(x, h):
    return (a(x,h)/(h*h))
def B_(x, h):
    return (a(x+h,h)/(h*h))
def C_(x, h):
    return((a(x+h,h)+a(x,h))/(h*h)+d(x,h))
def F_(x,h):
    return phi(x,h)
def K1_(h):
    return(a(h,h)/(h*(kappa0()+(h/2)*d0(h))+a(h,h)))
def V1_(h):
    return((h*(g0()+(h/2)*phi0(h)))/(h*(kappa0()+(h/2)*d0(h))+a(h,h)))
def Progonka_2(N):
    h = 1/N
    Y = [.0]*(N+1)
    alpha = [.0]*(N+1)
    beta = [.0]*(N+1)
    alpha[1] = K1_(h)
    beta[1] = V1_(h)
    for i in range(1, N):
        x_i = i*h
        alpha[i+1] = B_(x_i,h)/(C_(x_i, h)-A_(x_i, h)*alpha[i])
        beta[i+1] = (F_(x_i, h)+beta[i]*A_(x_i, h))/(C_(x_i, h)-A_(x_i, h)*alpha[i])
    Y[N] = g1()/kappa1()

    for i in range((N-1), -1, -1):
        Y[i] = alpha[i+1]*Y[i+1]+beta[i+1]
    return (Y)
#################3
def a_(x, h):
    return ((K(x)+K(x-h))/2)
def phi_(x,h):
    return (F(x))
def d_(x,h):
    return (q(x))
def d0_(h):
    return (q(0))
def phi0_(h):
    return (F(0))

def A__(x, h):
    return (a_(x,h)/(h*h))
def B__(x, h):
    return (a_(x+h,h)/(h*h))
def C__(x, h):
    return((a_(x+h,h)+a_(x,h))/(h*h)+d_(x,h))
def F__(x,h):
    return phi_(x,h)
def K1__(h):
    return(a_(h,h)/(h*(kappa0()+(h/2)*d0_(h))+a_(h,h)))
def V1__(h):
    return((h*(g0()+(h/2)*phi0_(h)))/(h*(kappa0()+(h/2)*d0_(h))+a_(h,h)))
def Progonka_3(N):
    h = 1/N
    Y = [.0]*(N+1)
    alpha = [.0]*(N+1)
    beta = [.0]*(N+1)
    alpha[1] = K1__(h)
    beta[1] = V1__(h)
    for i in range(1, N):
        x_i = i*h
        alpha[i+1] = B__(x_i,h)/(C__(x_i, h)-A__(x_i, h)*alpha[i])
        beta[i+1] = (F__(x_i, h)+beta[i]*A__(x_i, h))/(C__(x_i, h)-A__(x_i, h)*alpha[i])
    Y[N] = g1()/kappa1()

    for i in range((N-1), -1, -1):
        Y[i] = alpha[i+1]*Y[i+1]+beta[i+1]
    return (Y)
##################
def Print_Res(N1, N2):
    Y_t1000 = Progonka_1(N2)
    Y_t = [.0] * (N1+1)
    for i in range(N1+1):
        Y_t[i] = Y_t1000[i*(N2//N1)]
    print("Аппроксимация разностными операторами и повышение порядка точности:")
    Y1 = Progonka_1(N1)
    print("Значения Y1:")
    print_Y(Y1)
    print("Невязка точного решения с полученным приближенным:")
    print_pogr(Y_t, Y1)
    print("Интегро-интерполяционный метод:")
    Y2 = Progonka_2(N1)
    print("Значения Y2:")
    print_Y(Y2)
    print("Невязка точного решения с полученным приближенным:")
    print_pogr(Y_t, Y2)
    print("Вариационно-разностный метод:")
    Y3 = Progonka_3(N1)
    print("Значения Y3:")
    print_Y(Y3)
    print("Невязка точного решения с полученным приближенным:")
    print_pogr(Y_t, Y3)
###################main
Print_Res(10,1000)