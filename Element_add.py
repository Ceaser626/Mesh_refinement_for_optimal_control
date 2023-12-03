"""
This program adds the finite-element
"""

import numpy as np
import Parameters as para


def add_element(n_fe_list, n_cp_list, n_itv, d_list):
    """
    :param n_fe_list:
    :param n_cp_list:
    :param n_itv:
    :param d_list:
    :return:
    """

    """
    Initialize
    """
    count = 0

    for i in range(n_itv):
        """
        Get defection list
        """
        cp_sum = n_fe_list[i] * n_cp_list[i]
        d_itv_list = d_list[count:count + cp_sum]
        count = count + cp_sum

        """
        Calculate average error
        """
        d_aver = sum(d_itv_list) / len(d_itv_list)

        """
        Add collocation points
        """
        if d_aver >= para.tol_def:
            n_cp = int(np.ceil(n_cp_list[i] * (d_aver / para.tol_def) ** (1 / (n_cp_list[i] - 2.5))))
            if n_cp >= para.n_cp_max:
                n_fe = int(np.ceil(n_fe_list[i] * (d_aver / para.tol_def) ** (1 / (n_cp_list[i] - 2.5))))
                n_fe = min(n_fe, n_fe_list[i] + 3)  # At most add three
                n_fe_list[i] = n_fe
            else:
                n_cp_list[i] = n_cp

    return n_fe_list, n_cp_list, n_itv
