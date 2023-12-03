"""
This program adds the intervals and integrates the intervals
"""

import numpy as np
import Parameters as para


#################################################################################
#                        func: detect jumping intervals                         #
#################################################################################
def detect_jumping(ncol_list, d_list):
    """
    :param ncol_list:
    :param d_list:
    :return:
    """

    """
    Initialize
    """
    t_new_list = []
    l_new_list = []
    fix_flag_new_list = []

    """
    Average error
    """
    d_aver = sum(d_list) / len(d_list)
    d_aver_list = [d / d_aver for d in d_list]

    """
    Detect jumping interval
    """
    d_aver_max = 0
    for i in range(len(ncol_list)):
        if d_aver_list[i] >= para.tol_rel_def and d_list[i] >= para.tol_def:
            # Left point
            if d_aver_max == 0:
                t_left = ncol_list[i]
            # Maximum point
            if d_aver_list[i] >= d_aver_max:
                d_aver_max = d_aver_list[i]
                t_new = ncol_list[i]
        else:
            if d_aver_max != 0:
                t_right = ncol_list[i - 1]
                l_new = t_right - t_left
                t_new_list.append(t_new)
                l_new_list.append(l_new)
                fix_flag_new_list.append(0)
                d_aver_max = 0

    return t_new_list, l_new_list, fix_flag_new_list


#################################################################################
#                             func: divide intervals                            #
#################################################################################
def divide_interval(n_fe_list, n_cp_list, n_itv, output, t_fix_list, l_fix_list, fix_flag_list, t_new_list, l_new_list, fix_flag_new_list):
    """
    :param n_fe_list:
    :param n_cp_list:
    :param n_itv:
    :param output:
    :param t_fix_list:
    :param l_fix_list:
    :param fix_flag_list:
    :param t_new_list:
    :param l_new_list:
    :param fix_flag_new_list:
    :return:
    """

    """
    Initialize
    """
    # noncol_num = col_sum
    t_f_list = output[5]

    """
    Combine intervals
    """
    if t_fix_list != 0:
        t_add_list = t_new_list
        l_add_list = l_new_list
        fix_flag_add_list = fix_flag_new_list
        t_guess_list = t_fix_list + t_add_list
        l_guess_list = l_fix_list + l_add_list
        fix_flag_list = fix_flag_list + fix_flag_add_list
        # Rearrange the list
        zipped = zip(t_guess_list, l_guess_list, fix_flag_list)
        sort_zipped = sorted(zipped, key=lambda x: (x[0], x[1]))
        result = zip(*sort_zipped)
        t_guess_list, l_guess_list, fix_flag_list = [list(x) for x in result]
    else:
        t_add_list = t_new_list
        l_add_list = l_new_list
        t_guess_list = t_add_list
        l_guess_list = l_add_list
        fix_flag_list = fix_flag_new_list

    """
    Count new jump points in each interval
    """
    jump_num = []
    if t_fix_list != 0:
        for i in range(len(t_fix_list) + 1):
            if i == 0:
                t_fix_left = t_fix_list[i]
                num = sum(i < t_fix_left for i in t_guess_list)
            elif i == len(t_fix_list):
                t_fix_right = t_fix_list[i - 1]
                num = sum(i > t_fix_right for i in t_guess_list)
            else:
                t_fix_left = t_fix_list[i - 1]
                t_fix_right = t_fix_list[i]
                num = t_guess_list.index(t_fix_right) - t_guess_list.index(t_fix_left) - 1
            jump_num.append(num)
    else:
        jump_num.append(len(t_guess_list))

    """
    Divide interval
    """
    count = 0
    t_sum = 0
    n_fe_list_new = []
    n_cp_list_new = []
    for i in range(n_itv):
        col_sum = n_fe_list[i] * n_cp_list[i]
        n_cp_default = n_cp_list[i]
        itv_new = jump_num[i] + 1
        t_down = t_sum
        t_f = t_f_list[i]
        if itv_new > 1:
            for j in range(itv_new):
                if j != itv_new - 1:
                    ratio = (t_add_list[count] - t_down) / t_f
                    t_down = t_add_list[count]
                    count = count + 1
                else:
                    ratio = (t_sum + t_f - t_down) / t_f
                n_fe = int(np.ceil(ratio * col_sum / n_cp_default))
                n_cp = n_cp_default
                n_fe_list_new.append(n_fe)
                n_cp_list_new.append(n_cp)
        else:
            n_fe = n_fe_list[i]
            n_cp = n_cp_list[i]
            n_fe_list_new.append(n_fe)
            n_cp_list_new.append(n_cp)
        t_sum = t_sum + t_f

    n_itv_new = len(n_fe_list_new)

    return n_fe_list_new, n_cp_list_new, n_itv_new, t_guess_list, l_guess_list, fix_flag_list
