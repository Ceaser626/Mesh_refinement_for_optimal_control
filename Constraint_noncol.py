"""
This program computes the constraint at non-collocation points
"""

from Lagrange_interp import var_interp
import Parameters as para


def constraint(n_cp, data, tao):
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
    Compute constraint value
    """
    u_tao = U_tao[0]
    if 0 - 3 * para.tol_con <= u_tao <= 3 + 3 * para.tol_con:
        con_vil = 0
    else:
        con_vil = 1

    return con_vil
