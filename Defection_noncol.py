"""
This program computes the defection at non-collocation points
"""

import numpy as np
import Parameters as para
from Lagrange_interp import var_interp


def defection(n_cp, data, tao):
    """
    :param n_cp:
    :param data:
    :param tao:
    :return:
    """

    """
    Get interpolation value
    """
    X_tao, U_tao, X_dot_tao = var_interp(n_cp, data, tao)

    """
    Compute defection
    """
    h_tao = X_tao[0]
    v_tao = X_tao[1]
    m_tao = X_tao[2]
    u_tao = U_tao[0]
    fx = np.array([v_tao,
                   - para.g + u_tao,
                   - u_tao])
    t_f = data[4]

    defect = X_dot_tao - t_f * fx

    return defect
