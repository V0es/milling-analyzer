import time

import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

from numba.extending import overload, register_jitable
from numba import njit, complex128, void, float64

@register_jitable
def npexm(A):
    B = A.astype(np.complex128)
    eigvals, eigvecs = np.linalg.eig(B)
    diag_exp = np.diag(np.exp(eigvals))
    return (eigvecs @ diag_exp @ np.linalg.inv(eigvecs)).astype(complex128)

@overload(expm)
def numpy_expm(A):
    def expM(A):
        return npexm(A)
    return expM

# параметры
N = 2  # количество зубьев
Kt = 6e8  # коэффициент тангенциальной силы резания
Kn = 2e8  # коэффициент нормальной силы резания
w0 = 922 * 2 * np.pi  # собственная угловая частота (рад/с)
zeta = 0.011  # относительное затухание
m_t = 0.03993  # масса (кг)
aD = 0.05  # радиальная глубина резания
up_or_down = -1  # 1: встречное фрезерование, −1: попутное фрезерование

# углы входа и выхода
if up_or_down == 1:
    fist = 0  # угол входа
    fiex = np.arccos(1 - 2 * aD)  # угол выхода
elif up_or_down == -1:
    fist = np.arccos(2 * aD - 1)  # угол входа
    fiex = np.pi  # угол выхода

stx = 400  # шаги скорости шпинделя
sty = 200  # шаги глубины резания
w_st = 0e-3  # начальная глубина резания (м)
w_fi = 10e-3  # конечная глубина резания (м)
o_st = 5e3  # начальная скорость шпинделя (об/мин)
o_fi = 25e3  # конечная скорость шпинделя (об/мин)

# вычислительные параметры
k = 40  # количество интервалов дискретизации
intk = 20  # количество шагов численного интегрирования
m = k
wa = 0.5  # параметр времени задержки
wb = 0.5  # параметр времени задержки


# численное интегрирование для коэффициента резания
h_i = np.zeros(k)
for i in range(k):
    dtr = 2 * np.pi / N / k
    for j in range(1, N + 1):
        for h in range(intk):
            fi = i * dtr + (j - 1) * 2 * np.pi / N + h * dtr / intk
            g = 1 if fist <= fi <= fiex else 0
            h_i[i] += g * (Kt * np.cos(fi) + Kn * np.sin(fi)) * np.sin(fi) / intk


# начало выч
@njit()
def count():
    D = np.zeros((m + 2, m + 2), dtype=complex128)  # матрица D
    d = np.ones(m + 1, dtype=complex128)
    d[:2] = 0
    D += np.diag(d, -1)
    D[2, 0] = 1
    ss = np.zeros((stx + 1, sty + 1), dtype=complex128)
    dc = np.zeros((stx + 1, sty + 1), dtype=complex128)
    ei = np.zeros((stx + 1, sty + 1), dtype=np.complex128)
    for x in range(stx + 1):
        o = o_st + x * (o_fi - o_st) / stx  # скорость шпинделя
        tau = 60 / o / N  # время задержки
        dt = tau / m  # шаг времени

        for y in range(sty + 1):
            w = w_st + y * (w_fi - w_st) / sty  # глубина резания
            Fi = np.eye(m + 2, dtype=complex128)

            for i in range(m):
                A = np.zeros((2, 2), dtype=complex128)
                A[0, 1] = 1
                A[1, 0] = -w0 ** 2 - h_i[i] * w / m_t
                A[1, 1] = -2 * zeta * w0
                B = np.zeros((2, 2), dtype=complex128)
                B[1, 0] = h_i[i] * w / m_t

                P = npexm(A * dt)
                R = np.dot((npexm(A * dt) - np.eye(2, dtype=complex128)), np.linalg.inv(A)).dot(B)
                D[:2, :2] = P
                D[:2, m] = wa * R[:, 0]
                D[:2, m + 1] = wb * R[:, 0]
                Fi = np.dot(D, Fi)

            ss[x, y] = o  # скорость шпинделя
            dc[x, y] = w  # глубина резания
            ei[x, y] = max(np.abs(np.linalg.eigvals(Fi)))  # максимальное значение собственных чисел
    return ss, dc, ei
# график

t0 = time.time()
spin_speeds, depth_cut, eige = count()
print(f'time = {time.time() - t0}')
plt.figure()
plt.contour(spin_speeds, depth_cut, eige, levels=[1], colors='k')
plt.xlabel('Spindle Speed (rpm)')
plt.ylabel('Depth of Cut (m)')
plt.title('Stability Contour')
plt.show()
