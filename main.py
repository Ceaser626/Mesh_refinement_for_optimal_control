import Parameters as para
from Case_opt import opt, opt_adjust
from Result_anal import anal_result
from Interval_add import detect_jumping, divide_interval
from Element_add import add_element
from Plot import plot_result


#################################################################################
#                               0. Initialization                               #
#################################################################################
"""
0.1 Discretization information
"""
n_fe_list = [para.n_fe]
n_cp_list = [para.n_cp]
n_itv = len(n_fe_list)
"""
0.2 Algorithm initialization
"""
refine_flag = 0
iter_num = 0
t_fix_list = 0
l_fix_list = 0
fix_flag_list = 0
initial_value = 0


#################################################################################
#                         1. Attempt without refinement                         #
#################################################################################
"""
1.1 Optimization
"""
output = opt(n_fe_list, n_cp_list, n_itv, t_fix_list, iter_num, initial_value)
initial_value = output
"""
1.2 Analyse results
"""
col_list, ncol_list, d_list, h_list, c_list, u_0 = anal_result(output, n_fe_list, n_cp_list, n_itv, iter_num, save_name=1)
print('-'*20, f'Attempt', '-'*20, '\n', f'Finite-element: {n_fe_list}; Collocation points: {n_cp_list}', '\n')
"""
1.3 Judge necessity of refinement
"""
for i in range(len(ncol_list)):
    if d_list[i] > para.tol_def:
        refine_flag = 1


while refine_flag == 1 and iter_num <= para.iter_max:

    iter_num = iter_num + 1
    print('-'*20, f'Iteration {iter_num}', '-'*20)

    #################################################################################
    #                                2. Add interval                                #
    #################################################################################
    """
    2.1 Add new interval
    """
    t_new_list, l_new_list, fix_flag_new_list = detect_jumping(ncol_list, d_list)
    """
    2.2 Divide interval
    """
    n_fe_list, n_cp_list, n_itv, t_guess_list, l_guess_list, fix_flag_list = \
        divide_interval(n_fe_list, n_cp_list, n_itv, output, t_fix_list, l_fix_list, fix_flag_list, t_new_list, l_new_list, fix_flag_new_list)
    print(f'Divide interval', '\n', ' '*5, f'Guessed points: {t_guess_list}')
    """
    2.3 Adjust interval
    """
    for i in range(len(fix_flag_list)):
        fix_flag_list[i] = 0
    output, t_fix_list, l_fix_list = \
        opt_adjust(n_fe_list, n_cp_list, n_itv, t_guess_list, l_guess_list, iter_num, initial_value)
    initial_value = output
    print(f'Adjust interval', '\n', ' '*5, f'Interval points: {t_fix_list}', '\n', ' '*5, f'Finite-element: {n_fe_list}; Collocation points: {n_cp_list}')
    """
    2.4 Analyse results
    """
    col_list, ncol_list, d_list, h_list, c_list, u_0 = anal_result(output, n_fe_list, n_cp_list, n_itv, iter_num, save_name=0)

    #################################################################################
    #                             3. Add finite-element                             #
    #################################################################################
    """
    3.1 Add finite-element
    """
    n_fe_list, n_cp_list, n_itv = add_element(n_fe_list, n_cp_list, n_itv, d_list)
    print(f'Add finite-element', '\n', ' '*5, f'Finite-element: {n_fe_list}; Collocation points: {n_cp_list}', '\n')
    """
    3.2 Optimization
    """
    output = opt(n_fe_list, n_cp_list, n_itv, t_fix_list, iter_num, initial_value)
    initial_value = output
    """
    3.3 Analyse results
    """
    col_list, ncol_list, d_list, h_list, c_list, u_0 = anal_result(output, n_fe_list, n_cp_list, n_itv, iter_num, save_name=1)

    #################################################################################
    #                          4. Judge whether to end loop                         #
    #################################################################################
    """
    4.1 Judge necessity of refinement
    """
    refine_flag = 0
    for i in range(len(ncol_list)):
        if d_list[i] > para.tol_def:
            refine_flag = 1

print('-'*20, f'End', '-'*20, '\n', f'Total collocation points: {len(col_list)}')

plot_result(iter_num=4)
