"""
This program optimize the problem and return evaluations
"""

import numpy as np
from Lagrange_interp import var_interp
from Defection_noncol import defection
from Hamilton_col import hamilton
from Constraint_noncol import constraint


#################################################################################
#                             func: analyse result                              #
#################################################################################
def anal_result(output, n_fe_list, n_cp_list, n_itv, iter_num, save_name):
    """
    :param output:
    :param n_fe_list:
    :param n_cp_list:
    :param n_itv:
    :param iter_num:
    :param save_name:
    :return:
    """

    """
    0. Initialize 
    """
    # Load result from output
    tao_all = output[0]
    t_all = output[1]
    X_all = output[2]
    U_all = output[3]
    Dual_all = output[4]
    t_f_all = output[5]
    # Define list for storage
    col_list = []
    ncol_list = []
    d_list = []
    h_list = []
    c_list = []
    # Counter
    count = 0
    t_init = 0

    """
    1. Iterate overall intervals 
    """
    for itv in range(n_itv):  # iterate over intervals
        # Load numbers of finite-element, collocation point in current interval
        n_fe = n_fe_list[itv]
        n_cp = n_cp_list[itv]
        # Slice result for current interval
        itv_start = count
        itv_end = count + n_fe * n_cp + 1
        tao_itv = tao_all[itv_start:itv_end]
        t_itv = t_all[itv_start:itv_end]
        X_itv = X_all[:, itv_start:itv_end]
        U_itv = U_all[:, itv_start:itv_end]
        Dual_itv = Dual_all[:, itv_start:itv_end - 1]
        t_f_itv = t_f_all[itv]
        # Update
        tao_itv[0] = 0
        count = count + n_fe * n_cp
        t_init = t_itv[0]
        for i in range(n_fe):  # iterate over finite elements in the interval
            """
            2.0  Organize data in finite-element
            """
            # Slice result for current finite-element
            fe_start = i * n_cp
            fe_end = i * n_cp + n_cp + 1
            tao_list = tao_itv[fe_start:fe_end]
            t_list = t_itv[fe_start:fe_end]
            X_list = X_itv[:, fe_start:fe_end]
            U_list = U_itv[:, fe_start:fe_end]
            Dual_list = Dual_itv[:, fe_start:fe_end - 1]
            t_f = t_f_itv
            data = [tao_list, X_list, U_list, Dual_list, t_f]
            # Get collocation and non-collocation points
            non_col_tao_list = []
            for j in range(1, n_cp + 1):
                # For storage
                col = t_init + t_f_itv * tao_list[j]
                col_list.append(col)
                ncol = t_init + t_f_itv * (tao_list[j - 1] + tao_list[j]) / 2
                ncol_list.append(ncol)
                # For the following use
                non_col_tao = (tao_list[j - 1] + tao_list[j]) / 2
                non_col_tao_list.append(non_col_tao)

            """
            2.1 Calculate defection
            """
            for tao in non_col_tao_list:
                defect = defection(n_cp, data, tao)
                d_list.append(np.linalg.norm(defect))

            """
            2.2 Calculate hamilton
            """
            for j in range(n_cp):
                dual_list = Dual_list[:, j]
                tao = tao_list[j + 1]
                h = hamilton(n_fe, n_cp, data, tao, dual_list, j)
                h_list.append(h)

            """
            2.3 Calculate constraint
            """
            for tao in non_col_tao_list:
                con = constraint(n_cp, data, tao)
                c_list.append(con)

            """
            2.4 Calculate the control at 0
            """
            if itv == 0 and i == 0:
                tao = 0
                X_tao, U_tao, X_dot_tao = var_interp(n_cp, data, tao)
                u_0 = U_tao

    """
    3. Save information
    """
    if save_name == 1:
        np.savez(f'Figure/Err_iter{iter_num}_result', col_list=col_list,
                 ncol_list=ncol_list, d_list=d_list, h_list=h_list, c_list=c_list, u_0=u_0)
    else:
        np.savez(f'Figure/Err_mid{iter_num}_result', col_list=col_list,
                 ncol_list=ncol_list, d_list=d_list, h_list=h_list, c_list=c_list, u_0=u_0)

    """
    4. Return value
    """
    return col_list, ncol_list, d_list, h_list, c_list, u_0
