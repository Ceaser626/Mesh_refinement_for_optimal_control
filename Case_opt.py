"""
This file optimize Case 2 with the multi-resolution finite-element method
"""

import numpy as np
from pyomo.environ import *
from pyomo.dae import *
from Parameters import *
from Lagrange_interp import var_interp
from scipy.interpolate import interp1d


#################################################################################
#                            func: original problem                             #
#################################################################################
def opt(n_fe_list, n_cp_list, n_itv, t_fix_list, iter_num, initial_value):
    """
    :param n_fe_list:
    :param n_cp_list:
    :param n_itv:
    :param t_fix_list:
    :param iter_num:
    :param initial_value:
    :return:
    """

    m = ConcreteModel(name='opt_case2')
    m.itv = RangeSet(n_itv)  # from 1 to n_itv

    """
    For initialization use
    """
    if iter_num != 0:
        t = initial_value[1]
        X = initial_value[2]
        U = initial_value[3]
        X_f = interp1d(t, X, fill_value="extrapolate")
        U_f = interp1d(t, U, fill_value="extrapolate")
        t_list = []
        X_list = []
        U_list = []
        t_list.append(0)
        X_list.append(X_f(0))
        U_list.append(U_f(0))
        for i in t_fix_list:
            t_list.append(i)
            X_list.append(X_f(i))
            U_list.append(U_f(i))
        t_list.append(t[-1])
        X_list.append(X_f(t[-1]))
        U_list.append(U_f(t[-1]))

    """
    Create a block for a single interval
    """
    def interval_block_rule(c, itv):
        ##############################
        #      Define variables      #
        ##############################
        n_fe = n_fe_list[itv - 1]
        n_cp = n_cp_list[itv - 1]
        # Time
        c.tau = ContinuousSet(bounds=(0, 1))
        c.t_f = Var(bounds=(0, 5))
        c.t = Var(c.tau, bounds=(0, 5))
        # State variable
        c.h = Var(c.tau, bounds=(0, 10))
        c.v = Var(c.tau, bounds=(-5, 0))
        c.mass = Var(c.tau, bounds=(0, 10))
        c.u = Var(c.tau, bounds=(0, 3), initialize=1)
        # Corresponding state variable
        c.d_h = DerivativeVar(c.h)
        c.d_v = DerivativeVar(c.v)
        c.d_mass = DerivativeVar(c.mass)
        c.d_t = DerivativeVar(c.t)

        ##############################
        #       Discretization       #
        ##############################
        discretizer = TransformationFactory('dae.collocation')
        discretizer.apply_to(c, wrt=c.tau, nfe=n_fe, ncp=n_cp, scheme='LAGRANGE-RADAU')

        ##############################
        #       Initialization       #
        ##############################
        if iter_num != 0:
            X_init = X_list[itv - 1]
            X_fina = X_list[itv]
            U_init = U_list[itv - 1]
            U_fina = U_list[itv]
            c.t_f = np.clip(t_list[itv] - t_list[itv - 1], 0, 5)
            for i in c.tau:
                c.t[i] = np.clip((1 - value(i)) * t_list[itv - 1] + value(i) * t_list[itv], 0, 5)
                c.h[i] = np.clip((1 - value(i)) * X_init[0] + value(i) * X_fina[0], 0, 10)
                c.v[i] = np.clip((1 - value(i)) * X_init[1] + value(i) * X_fina[1], -5, 0)
                c.mass[i] = np.clip((1 - value(i)) * X_init[2] + value(i) * X_fina[2], 0, 10)
                c.u[i] = np.clip((1 - value(i)) * U_init[0] + value(i) * U_fina[0], 0, 3)

        ##############################
        # Define dynamic constraints #
        ##############################
        def d_h_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_h[i] == c.t_f * c.v[i]
        c.ode_h = Constraint(c.tau, rule=d_h_rule)

        def d_v_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_v[i] == c.t_f * (- g + c.u[i])
        c.ode_v = Constraint(c.tau, rule=d_v_rule)

        def d_mass_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_mass[i] == c.t_f * (- c.u[i])
        c.ode_mass = Constraint(c.tau, rule=d_mass_rule)

        def d_t_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_t[i] == c.t_f
        c.ode_t = Constraint(c.tau, rule=d_t_rule)

        ##############################
        #   Define path constraints  #
        ##############################
        def pathcon_u_low_rule(c, i):
            return 0 <= c.u[i]
        c.pathcon_u_low = Constraint(c.tau, rule=pathcon_u_low_rule)

        def pathcon_u_up_rule(c, i):
            return c.u[i] <= 3
        c.pathcon_u_up = Constraint(c.tau, rule=pathcon_u_up_rule)

    m.inst = Block(m.itv, rule=interval_block_rule)

    """
    Define boundary constraints
    """
    def bound_h_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].h[0] == 10
        else:
            return Constraint.Skip
    m.bound_h_i = Constraint(m.itv, rule=bound_h_i_rule)

    def bound_v_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].v[0] == - 2
        else:
            return Constraint.Skip
    m.bound_v_i = Constraint(m.itv, rule=bound_v_i_rule)

    def bound_mass_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].mass[0] == 10
        else:
            return Constraint.Skip
    m.bound_mass_i = Constraint(m.itv, rule=bound_mass_i_rule)

    def bound_t_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].t[0] == 0
        else:
            return Constraint.Skip
    m.bound_t_i = Constraint(m.itv, rule=bound_t_i_rule)

    def bound_h_f_rule(m, itv):
        if itv == m.itv.last():
            return m.inst[itv].h[1] == 0
        else:
            return Constraint.Skip
    m.bound_h_f = Constraint(m.itv, rule=bound_h_f_rule)

    def bound_v_f_rule(m, itv):
        if itv == m.itv.last():
            return m.inst[itv].v[1] == 0
        else:
            return Constraint.Skip
    m.bound_v_f = Constraint(m.itv, rule=bound_v_f_rule)

    """
    Define linking constraints
    """
    def linking_h_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].h[0] == m.inst[itv - 1].h[1]
        else:
            return Constraint.Skip
    m.linking_h = Constraint(m.itv, rule=linking_h_rule)

    def linking_v_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].v[0] == m.inst[itv - 1].v[1]
        else:
            return Constraint.Skip
    m.linking_v = Constraint(m.itv, rule=linking_v_rule)

    def linking_mass_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].mass[0] == m.inst[itv - 1].mass[1]
        else:
            return Constraint.Skip
    m.linking_mass = Constraint(m.itv, rule=linking_mass_rule)

    def linking_t_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].t[0] == m.inst[itv - 1].t[1]
        else:
            return Constraint.Skip
    m.linking_t = Constraint(m.itv, rule=linking_t_rule)

    """
    Define jump position constraints
    """
    def jump_t_rule(m, itv):
        if itv != m.itv.last() and t_fix_list != 0:
            return m.inst[itv].t[1] == t_fix_list[itv - 1]
        else:
            return Constraint.Skip
    m.jump_t = Constraint(m.itv, rule=jump_t_rule)

    """
    Define cost function
    """
    def objective_rule(m):
        itv = m.itv.last()
        return - m.inst[itv].mass[1]
    m.objective = Objective(rule=objective_rule, sense=minimize)

    """
    Solve
    """
    solver = SolverFactory('ipopt', sol='nl')
    solver.reset()
    solver.options['halt_on_ampl_error'] = "yes"
    solver.options['linear_solver'] = "ma57"
    solver.options['tol'] = "1e-6"
    solver.options['mu_strategy'] = "adaptive"
    solver.options['bound_push'] = "1e-6"
    solver.options['max_iter'] = 50000
    m.dual = Suffix(direction=Suffix.IMPORT)
    try:
        results = solver.solve(m, tee=False)
        term_cond = results.solver.termination_condition
        if term_cond == TerminationCondition.optimal:
            cond = 'Success'
        else:
            cond = 'Not Success'
    except Exception:
        cond = 'Failure'

    """
    Store
    """
    tao = []
    t = []
    h = []
    v = []
    mass = []
    u = []
    t_f = []
    dual_ode_h = []
    dual_ode_v = []
    dual_ode_m = []
    for i in m.inst:
        for j in m.inst[i].tau:
            if i == 1 or j != 0:
                tao.append(j)
                t.append(value(m.inst[i].t[j]))
                h.append(value(m.inst[i].h[j]))
                v.append(value(m.inst[i].v[j]))
                mass.append(value(m.inst[i].mass[j]))
                u.append(value(m.inst[i].u[j]))
        t_f.append(value(m.inst[i].t_f))
    for i in m.inst:
        for index in m.inst[i].ode_h:
            dual_ode_h.append(m.dual[m.inst[i].ode_h[index]])
            dual_ode_v.append(m.dual[m.inst[i].ode_v[index]])
            dual_ode_m.append(m.dual[m.inst[i].ode_mass[index]])
    np.savez(f'Figure/Sol_iter{iter_num}_result',
             tao=tao, t=t, h=h, v=v, mass=mass, u=u, t_f=t_f, n_fe_list=n_fe_list, n_cp_list=n_cp_list)

    """
    Output
    """
    tao = np.array(tao)
    t = np.array(t)
    h = np.array(h)
    v = np.array(v)
    mass = np.array(mass)
    u = np.array(u)
    dual_ode_h = np.array(dual_ode_h)
    dual_ode_v = np.array(dual_ode_v)
    dual_ode_m = np.array(dual_ode_m)
    X = np.vstack([h, v, mass])
    U = np.vstack([u])
    Dual = np.vstack([dual_ode_h, dual_ode_v, dual_ode_m])
    output = [tao, t, X, U, Dual, t_f]

    return output


#################################################################################
#                         func: adjust interval length                          #
#################################################################################
def opt_adjust(n_fe_list, n_cp_list, n_itv, t_guess_list, l_guess_list, iter_num, initial_value):
    """
    :param n_fe_list:
    :param n_cp_list:
    :param n_itv:
    :param t_guess_list:
    :param l_guess_list:
    :param iter_num:
    :param initial_value:
    :return:
    """

    m = ConcreteModel(name='opt_case2_adjust')
    m.itv = RangeSet(n_itv)  # from 1 to n_itv

    """
    For initialization use
    """
    t = initial_value[1]
    X = initial_value[2]
    U = initial_value[3]
    X_f = interp1d(t, X, fill_value="extrapolate")
    U_f = interp1d(t, U, fill_value="extrapolate")
    t_list = []
    X_list = []
    U_list = []
    t_list.append(0)
    X_list.append(X_f(0))
    U_list.append(U_f(0))
    for i in t_guess_list:
        t_list.append(i)
        X_list.append(X_f(i))
        U_list.append(U_f(i))
    t_list.append(t[-1])
    X_list.append(X_f(t[-1]))
    U_list.append(U_f(t[-1]))

    """
    Create a block for a single interval
    """
    def interval_block_rule(c, itv):
        ##############################
        #      Define variables      #
        ##############################
        n_fe = n_fe_list[itv - 1]
        n_cp = n_cp_list[itv - 1]
        n_var = 3
        c.n_fe = RangeSet(n_fe)
        c.n_cp = RangeSet(n_cp)
        c.n_u_cp = RangeSet(n_cp - 1)
        c.n_var = RangeSet(n_var)
        # Time
        c.tau = ContinuousSet(bounds=(0, 1))
        c.t_f = Var(bounds=(0.1, 5))
        c.t = Var(c.tau, bounds=(0, 5))
        # State variable
        c.h = Var(c.tau, bounds=(0, 10))
        c.v = Var(c.tau, bounds=(-5, 0))
        c.mass = Var(c.tau, bounds=(0, 10))
        c.u = Var(c.tau, bounds=(0, 3), initialize=1)
        # Corresponding state variable
        c.d_h = DerivativeVar(c.h)
        c.d_v = DerivativeVar(c.v)
        c.d_mass = DerivativeVar(c.mass)
        c.d_t = DerivativeVar(c.t)
        # Error var
        c.err = Var(c.n_fe, c.n_cp, c.n_var)
        c.err_sum = Var()

        ##############################
        #       Discretization       #
        ##############################
        discretizer = TransformationFactory('dae.collocation')
        discretizer.apply_to(c, wrt=c.tau, nfe=n_fe, ncp=n_cp, scheme='LAGRANGE-RADAU')

        ##############################
        #       Initialization       #
        ##############################
        X_init = X_list[itv - 1]
        X_fina = X_list[itv]
        U_init = U_list[itv - 1]
        U_fina = U_list[itv]
        c.t_f = np.clip(t_list[itv] - t_list[itv - 1], 0, 5)
        for i in c.tau:
            c.t[i] = np.clip((1 - value(i)) * t_list[itv - 1] + value(i) * t_list[itv], 0, 5)
            c.h[i] = np.clip((1 - value(i)) * X_init[0] + value(i) * X_fina[0], 0, 10)
            c.v[i] = np.clip((1 - value(i)) * X_init[1] + value(i) * X_fina[1], -5, 0)
            c.mass[i] = np.clip((1 - value(i)) * X_init[2] + value(i) * X_fina[2], 0, 10)
            c.u[i] = np.clip((1 - value(i)) * U_init[0] + value(i) * U_fina[0], 0, 3)

        ##############################
        # Define dynamic constraints #
        ##############################
        def d_h_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_h[i] == c.t_f * c.v[i]
        c.ode_h = Constraint(c.tau, rule=d_h_rule)

        def d_v_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_v[i] == c.t_f * (- g + c.u[i])
        c.ode_v = Constraint(c.tau, rule=d_v_rule)

        def d_mass_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_mass[i] == c.t_f * (- c.u[i])
        c.ode_mass = Constraint(c.tau, rule=d_mass_rule)

        def d_t_rule(c, i):
            if i == 0:
                return Constraint.Skip
            else:
                return c.d_t[i] == c.t_f
        c.ode_t = Constraint(c.tau, rule=d_t_rule)

        #####################################
        #  Define control sign constraints  #
        #####################################
        def pathcon_control_sign_rule(c, i, j):
            # Get var in finite-element
            tao_list = []
            u_list = []
            for l in range(1, n_cp + 1):
                tao_l = c.tau.at((i - 1) * n_cp + l + 1)
                tao_list.append(tao_l)
                u_list.append(c.u[tao_l])
            tao_cur = c.tau.at((i - 1) * n_cp + j + 1)
            tao_nex = c.tau.at((i - 1) * n_cp + j + 2)
            # Interpolate at non-collocation points
            X_list = np.vstack([u_list])
            U_list = np.vstack([u_list])
            data = [tao_list, X_list, U_list]
            X_tao, U_tao, X_dot_tao_cur = var_interp(n_cp - 1, data, tao_cur)
            X_tao, U_tao, X_dot_tao_nex = var_interp(n_cp - 1, data, tao_nex)
            u_dot_tao_cur = X_dot_tao_cur[0]
            u_dot_tao_nex = X_dot_tao_nex[0]
            return - u_dot_tao_cur * u_dot_tao_nex <= 0
        c.pathcon_control_sign = Constraint(c.n_fe, c.n_u_cp, rule=pathcon_control_sign_rule)

        ##############################
        #  Define error constraints  #
        ##############################
        def bound_err_rule(c, i, j, k):
            return c.err[i, j, k] >= 0
        c.bound_err = Constraint(c.n_fe, c.n_cp, c.n_var, rule=bound_err_rule)

        def pathcon_err_low_rule(c, i, j, k):
            # Get var in finite-element
            tao_list = []
            h_list = []
            v_list = []
            mass_list = []
            u_list = []
            for l in range(n_cp + 1):
                tao_l = c.tau.at((i - 1) * n_cp + l + 1)
                tao_list.append(tao_l)
                h_list.append(c.h[tao_l])
                v_list.append(c.v[tao_l])
                mass_list.append(c.mass[tao_l])
                u_list.append(c.u[tao_l])
            tao = (c.tau.at((i - 1) * n_cp + j) + c.tau.at((i - 1) * n_cp + j + 1)) / 2
            # Interpolate at non-collocation points
            X_list = np.vstack([h_list, v_list, mass_list])
            U_list = np.vstack([u_list])
            data = [tao_list, X_list, U_list]
            X_tao, U_tao, X_dot_tao = var_interp(n_cp, data, tao)
            h_tao = X_tao[0]
            v_tao = X_tao[1]
            mass_tao = X_tao[2]
            u_tao = U_tao[0]
            h_dot_tao = X_dot_tao[0]
            v_dot_tao = X_dot_tao[1]
            mass_dot_tao = X_dot_tao[2]
            # Get error
            if k == 1:
                err = h_dot_tao - c.t_f * v_tao
            elif k == 2:
                err = v_dot_tao - c.t_f * (- g + u_tao)
            elif k == 3:
                err = mass_dot_tao - c.t_f * (- u_tao)
            return - c.err[i, j, k] <= err
        c.pathcon_err_low = Constraint(c.n_fe, c.n_cp, c.n_var, rule=pathcon_err_low_rule)

        def pathcon_err_up_rule(c, i, j, k):
            # Get var in finite-element
            tao_list = []
            h_list = []
            v_list = []
            mass_list = []
            u_list = []
            for l in range(n_cp + 1):
                tao_l = c.tau.at((i - 1) * n_cp + l + 1)
                tao_list.append(tao_l)
                h_list.append(c.h[tao_l])
                v_list.append(c.v[tao_l])
                mass_list.append(c.mass[tao_l])
                u_list.append(c.u[tao_l])
            tao = (c.tau.at((i - 1) * n_cp + j) + c.tau.at((i - 1) * n_cp + j + 1)) / 2
            # Interpolate at non-collocation points
            X_list = np.vstack([h_list, v_list, mass_list])
            U_list = np.vstack([u_list])
            data = [tao_list, X_list, U_list]
            X_tao, U_tao, X_dot_tao = var_interp(n_cp, data, tao)
            h_tao = X_tao[0]
            v_tao = X_tao[1]
            mass_tao = X_tao[2]
            u_tao = U_tao[0]
            h_dot_tao = X_dot_tao[0]
            v_dot_tao = X_dot_tao[1]
            mass_dot_tao = X_dot_tao[2]
            # Get error
            if k == 1:
                err = h_dot_tao - c.t_f * v_tao
            elif k == 2:
                err = v_dot_tao - c.t_f * (- g + u_tao)
            elif k == 3:
                err = mass_dot_tao - c.t_f * (- u_tao)
            return err <= c.err[i, j, k]
        c.pathcon_con_err_up = Constraint(c.n_fe, c.n_cp, c.n_var, rule=pathcon_err_up_rule)

        def pathcon_err_sum_rule(c):
            return c.err_sum == sum(c.err[i, j, k] for i in c.n_fe for j in c.n_cp for k in c.n_var)
        c.pathcon_err_sum = Constraint(rule=pathcon_err_sum_rule)

        ##############################
        #   Define path constraints  #
        ##############################
        def pathcon_u_low_rule(c, i):
            return 0 <= c.u[i]
        c.pathcon_u_low = Constraint(c.tau, rule=pathcon_u_low_rule)

        def pathcon_u_up_rule(c, i):
            return c.u[i] <= 3
        c.pathcon_u_up = Constraint(c.tau, rule=pathcon_u_up_rule)

    m.inst = Block(m.itv, rule=interval_block_rule)

    """
    Define boundary constraints
    """

    def bound_h_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].h[0] == 10
        else:
            return Constraint.Skip
    m.bound_h_i = Constraint(m.itv, rule=bound_h_i_rule)

    def bound_v_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].v[0] == - 2
        else:
            return Constraint.Skip
    m.bound_v_i = Constraint(m.itv, rule=bound_v_i_rule)

    def bound_mass_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].mass[0] == 10
        else:
            return Constraint.Skip
    m.bound_mass_i = Constraint(m.itv, rule=bound_mass_i_rule)

    def bound_t_i_rule(m, itv):
        if itv == m.itv.first():
            return m.inst[itv].t[0] == 0
        else:
            return Constraint.Skip
    m.bound_t_i = Constraint(m.itv, rule=bound_t_i_rule)

    def bound_h_f_rule(m, itv):
        if itv == m.itv.last():
            return m.inst[itv].h[1] == 0
        else:
            return Constraint.Skip
    m.bound_h_f = Constraint(m.itv, rule=bound_h_f_rule)

    def bound_v_f_rule(m, itv):
        if itv == m.itv.last():
            return m.inst[itv].v[1] == 0
        else:
            return Constraint.Skip
    m.bound_v_f = Constraint(m.itv, rule=bound_v_f_rule)

    """
    Define linking constraints
    """
    def linking_h_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].h[0] == m.inst[itv - 1].h[1]
        else:
            return Constraint.Skip
    m.linking_h = Constraint(m.itv, rule=linking_h_rule)

    def linking_v_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].v[0] == m.inst[itv - 1].v[1]
        else:
            return Constraint.Skip
    m.linking_v = Constraint(m.itv, rule=linking_v_rule)

    def linking_mass_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].mass[0] == m.inst[itv - 1].mass[1]
        else:
            return Constraint.Skip
    m.linking_mass = Constraint(m.itv, rule=linking_mass_rule)

    def linking_t_rule(m, itv):
        if itv != m.itv.first():
            return m.inst[itv].t[0] == m.inst[itv - 1].t[1]
        else:
            return Constraint.Skip
    m.linking_t = Constraint(m.itv, rule=linking_t_rule)

    """
    Define jump position constraints
    """
    def jump_t_rule(m, itv):
        if itv != m.itv.last():
            left_bound = (t_guess_list[itv - 1] - l_guess_list[itv - 1])
            right_bound = (t_guess_list[itv - 1] + l_guess_list[itv - 1])
            return inequality(left_bound, m.inst[itv].t[1], right_bound)
        else:
            return Constraint.Skip
    m.jump_t = Constraint(m.itv, rule=jump_t_rule)

    """
    Define cost function
    """
    def objective_rule(m):
        fina = m.itv.last()
        return - m.inst[fina].mass[1] + omiga_1 * sum(m.inst[itv].err_sum for itv in m.itv)
    m.objective = Objective(rule=objective_rule, sense=minimize)

    """
    Solve
    """
    solver = SolverFactory('ipopt', sol='nl')
    solver.reset()
    solver.options['halt_on_ampl_error'] = "yes"
    solver.options['linear_solver'] = "ma57"
    solver.options['tol'] = "1e-6"
    solver.options['mu_strategy'] = "adaptive"
    solver.options['bound_push'] = "1e-6"
    solver.options['max_iter'] = 50000
    m.dual = Suffix(direction=Suffix.IMPORT)
    try:
        results = solver.solve(m, tee=False)
        term_cond = results.solver.termination_condition
        if term_cond == TerminationCondition.optimal:
            cond = 'Success'
        else:
            cond = 'Not Success'
    except Exception:
        cond = 'Failure'

    """
    Store
    """
    tao = []
    t = []
    h = []
    v = []
    mass = []
    u = []
    t_f = []
    dual_ode_h = []
    dual_ode_v = []
    dual_ode_m = []
    t_fix_list = []
    l_fix_list = []
    count = 0
    for i in m.inst:
        for j in m.inst[i].tau:
            if i == 1 or j != 0:
                tao.append(j)
                t.append(value(m.inst[i].t[j]))
                h.append(value(m.inst[i].h[j]))
                v.append(value(m.inst[i].v[j]))
                mass.append(value(m.inst[i].mass[j]))
                u.append(value(m.inst[i].u[j]))
        t_f.append(value(m.inst[i].t_f))
        if i != m.itv.last():
            t_fix_list.append(value(m.inst[i].t[1]))
            l_fix_list.append(np.abs(value(m.inst[i].t[1]) - t_guess_list[count]))
            count = count + 1
    for i in m.inst:
        for index in m.inst[i].ode_h:
            dual_ode_h.append(m.dual[m.inst[i].ode_h[index]])
            dual_ode_v.append(m.dual[m.inst[i].ode_v[index]])
            dual_ode_m.append(m.dual[m.inst[i].ode_mass[index]])
    np.savez(f'Figure/Adjust_iter{iter_num}_result',
             tao=tao, t=t, h=h, v=v, mass=mass, u=u, t_f=t_f, n_fe_list=n_fe_list, n_cp_list=n_cp_list)

    """
    Output
    """
    tao = np.array(tao)
    t = np.array(t)
    h = np.array(h)
    v = np.array(v)
    mass = np.array(mass)
    u = np.array(u)
    dual_ode_h = np.array(dual_ode_h)
    dual_ode_v = np.array(dual_ode_v)
    dual_ode_m = np.array(dual_ode_m)
    X = np.vstack([h, v, mass])
    U = np.vstack([u])
    Dual = np.vstack([dual_ode_h, dual_ode_v, dual_ode_m])
    output = [tao, t, X, U, Dual, t_f]

    return output, t_fix_list, l_fix_list
