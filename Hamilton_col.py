"""
This program computes the Hamilton at collocation points
"""

import numpy as np
from Lagrange_interp import var_interp


def collocation_weight(n_cp):
    if n_cp == 2:
        weight = [1.5, 0.5]
    elif n_cp == 3:
        weight = [0.752806125400934, 1.024971652376843, 0.222222222222222]
    elif n_cp == 4:
        weight = [0.440924422353536, 0.776386937686344, 0.657688639960120, 0.125000000000000]
    elif n_cp == 5:
        weight = [0.287427121582451, 0.562712030298924, 0.623653045951482, 0.446207802167142, 0.080000000000000]
    elif n_cp == 6:
        weight = [0.201588385253480, 0.416901334311909, 0.520926783189575, 0.485387188468970, 0.319640753220511, 0.055555555555556]
    elif n_cp == 7:
        weight = [0.148988471112020, 0.318204231467302, 0.424703779005956, 0.447109829014567, 0.380949873644231, 0.239227489225312, 0.040816326530612]
    elif n_cp == 8:
        weight = [0.114508814744257, 0.249647901329864, 0.347014795634501, 0.391572167452494, 0.376517545389119, 0.304130620646785, 0.185358154802979, 0.031250000000000]
    elif n_cp == 9:
        weight = [0.090714504923287, 0.200553298024551, 0.286386696357231, 0.337693966975929, 0.348273002772967, 0.316843775670438, 0.247189378204593, 0.147654019046315, 0.024691358024691]
    elif n_cp == 10:
        weight = [0.073617005486761, 0.164376012736922, 0.239193431714380, 0.290610164832918, 0.313582457226938, 0.305859287724423, 0.268194837841179, 0.204270131879001, 0.120296670557482, 0.020000000000000]
    else:
        weight = np.zeros(n_cp)
    return weight


def hamilton(n_fe, n_cp, data, tao, dual_list, j):
    """
    :param n_cp:
    :param data:
    :param tao:
    :param dual_list:
    :param t_f:
    :param j:
    :param n_fe:
    :return:
    """

    """
    Get interpolation value
    """
    X_tao, U_tao, X_dot_tao = var_interp(n_cp, data, tao)

    """
    Get weight
    """
    weight = collocation_weight(n_cp)

    """
    Compute hamilton
    """
    t_f = data[4]
    h = np.inner(dual_list, X_dot_tao) / weight[j] / t_f * n_fe

    return h
