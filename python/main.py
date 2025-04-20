import time
import math
import logging
import numpy as np
import argparse as ap
import multiprocessing as mp
from itertools import repeat

logger = logging.getLogger(__name__)


SIZE_T = 1000
SIZE_X = 500

size = SIZE_X - 1
m = 2
T = 5
h = 0.02
g = 100
t = 0.01
r = t / h
N = 5


def cov_alpha(m: float, T: float, h: float):
    return (m * m / T) * h + 2 / (T * h)


def cov_beta(T: float, h: float):
    return -1 / T / h


def matrix_a():
    return np.diag(np.full(size, h))


def matrix_b():
    alpha = cov_alpha(m, T, h)
    beta = cov_beta(T, h)

    matrix = np.zeros((size, size))
    for i in range(size):
        matrix[i, i] = alpha
        if i + 1 < size:
            matrix[i, i + 1] = beta
        if i - 1 >= 0:
            matrix[i, i - 1] = beta

    matrix[0, size - 1] = beta
    matrix[size - 1, 0] = beta

    return matrix


def F(j: int, i: int, phi):
    return m * m * phi[j, i] * phi[j, i] / 2 + g * math.pow(phi[j, i], 4) / 4


def Strauss_Vazquez(phi):
    def boundary_condition_1():
        if phi[j + 1, i] == phi[j - 1, i]:
            phi[j + 1, i] = r * r * phi[j, SIZE_X - 2] + 2 * (1 - r * r) * phi[j, i] + r * r * phi[j, i + 1] - phi[j - 1, i] - t * t * (m * m * phi[j + 1, i] + g * math.pow(phi[j + 1, i], 3))
        else:
            phi[j + 1, i] = r * r * phi[j, SIZE_X - 2] + 2 * (1 - r * r) * phi[j, i] + r * r * phi[j, i + 1] - phi[j - 1, i] - t * t * (F(j + 1, i, phi) - F(j - 1, i, phi)) / (phi[j + 1, i] - phi[j - 1, i])

    def boundary_condition_2():
        if phi[j + 1, i] == phi[j - 1, i]:
            phi[j + 1, i] = r * r * phi[j, i - 1] + 2 * (1 - r * r) * phi[j, i] + r * r * phi[j, 1] - phi[j - 1, i] - t * t * (m * m * phi[j + 1, i] + g * math.pow(phi[j + 1, i], 3))
        else:
            phi[j + 1, i] = r * r * phi[j, i - 1] + 2 * (1 - r * r) * phi[j, i] + r * r * phi[j, 1] - phi[j - 1, i] - t * t * (F(j + 1, i, phi) - F(j - 1, i, phi)) / (phi[j + 1, i] - phi[j - 1, i])

    def boundary_condition_3():
        if phi[j + 1, i] == phi[j - 1, i]:
            phi[j + 1, i] = r * r * phi[j, i - 1] + 2 * (1 - r * r) * phi[j, i] + r * r * phi[j, i + 1] - phi[j - 1, i] - t * t * (m * m * phi[j + 1, i] + g * math.pow(phi[j + 1, i], 3))
        else:
            phi[j + 1, i] = r * r * phi[j, i - 1] + 2 * (1 - r * r) * phi[j, i] + r * r * phi[j, i + 1] - phi[j - 1, i] - t * t * (F(j + 1, i, phi) - F(j - 1, i, phi)) / (phi[j + 1, i] - phi[j - 1, i])

    for j in range(1, SIZE_T - 1):
        for i in range(0, SIZE_X):
            phi[j + 1, i] = phi[j, i]
            recent = phi[j + 1, i]

            if i == 0:
                cnt = 0
                while abs(recent - phi[j + 1, i]) > 1e-5 or cnt == 0:
                    recent = phi[j + 1, i]
                    boundary_condition_1()
                    cnt += 1
            elif i == SIZE_X - 1:
                cnt = 0
                while abs(recent - phi[j + 1, i]) > 1e-5 or cnt == 0:
                    recent = phi[j + 1, i]
                    boundary_condition_2()
                    cnt += 1
            else:
                cnt = 0
                while abs(recent - phi[j + 1, i]) > 1e-5 or cnt == 0:
                    recent = phi[j + 1, i]
                    boundary_condition_3()
                    cnt += 1


def main(mean, cov_a, cov_b):
    multivariate_normal_a = np.random.multivariate_normal(mean, cov_a)
    multivariate_normal_b = np.random.multivariate_normal(mean, cov_b)

    phi = np.zeros((SIZE_T, SIZE_X))
    phi[0, :-1] = multivariate_normal_b
    phi[1, :-1] = phi[0, :-1] + t * multivariate_normal_a
    phi[:, -1] = phi[:, 0]

    Strauss_Vazquez(phi)

    return phi


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    parser = ap.ArgumentParser()
    parser.add_argument("--jobs", "-j", type=int, default=1)

    args = parser.parse_args()
    jobs = args.jobs if args.jobs <= mp.cpu_count() else mp.cpu_count()
    logger.info("CPU used: %s" % jobs)

    mean = np.zeros(size, dtype=np.float64)
    cov_a = np.linalg.inv(matrix_a())
    cov_b = np.linalg.inv(matrix_b())

    mean_phi = np.zeros((SIZE_T, SIZE_X))

    # for _ in tqdm(range(N)):
    #     multivariate_normal_a = np.random.multivariate_normal(mean, cov_a)
    #     multivariate_normal_b = np.random.multivariate_normal(mean, cov_b)

    #     phi = np.zeros((SIZE_T, SIZE_X))
    #     phi[0, :-1] = multivariate_normal_b
    #     phi[1, :-1] = phi[0, :-1] + t * multivariate_normal_a
    #     phi[:, -1] = phi[:, 0]

    #     Strauss_Vazquez(phi)

    #     mean_phi += phi / N

    t_start = time.perf_counter()

    with mp.Pool(jobs) as pool:
        pool.starmap(main, zip(
            repeat(mean, N),
            repeat(cov_a, N),
            repeat(cov_b, N)
        ))

    t_end = time.perf_counter()
    logger.info("It took %s seconds" % str(t_end - t_start))
