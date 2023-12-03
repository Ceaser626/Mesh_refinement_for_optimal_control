"""
This program interpolate the state and control variable
"""
import numpy as np


#################################################################################
#            func: interpolate the variable from optimization result            #
#################################################################################
def var_interp(n_cp, data, tao):
    """
    :param n_cp:
    :param data:
    :param tao:
    :return:
    """

    """
    Get data
    """
    tao_list = data[0]
    X_list = data[1]
    U_list = data[2]

    """
    Interpolate X
    """
    X_tao = np.zeros(np.size(X_list[:, 0]))
    for j in range(n_cp + 1):
        cof = 1
        for k in range(n_cp + 1):
            if k != j:
                cof = cof * (tao - tao_list[k]) / (tao_list[j] - tao_list[k])
        X_tao = X_tao + cof * X_list[:, j]

    """
    Interpolate U
    """
    U_tao = np.zeros(np.size(U_list[:, 0]))
    for j in range(1, n_cp + 1):
        cof = 1
        for k in range(1, n_cp + 1):
            if k != j:
                cof = cof * (tao - tao_list[k]) / (tao_list[j] - tao_list[k])
        U_tao = U_tao + cof * U_list[:, j]

    """
    Interpolate X_dot
    """
    X_dot_tao = np.zeros(np.size(X_list[:, 0]))
    for j in range(n_cp + 1):
        cof_n = 0
        cof_d = 1
        cof = 1
        for s in range(n_cp + 1):
            if s != j:
                cof_in = 1
                for k in range(n_cp + 1):
                    if k != j and k != s:
                        cof_in = cof_in * (tao - tao_list[k])
                cof_n = cof_n + cof_in
        for k in range(n_cp + 1):
            if k != j:
                cof_d = cof_d * (tao_list[j] - tao_list[k])
        cof = cof * cof_n / cof_d
        X_dot_tao = X_dot_tao + cof * X_list[:, j]

    return X_tao, U_tao, X_dot_tao
