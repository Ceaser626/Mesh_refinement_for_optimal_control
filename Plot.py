import numpy as np
import matplotlib.pyplot as plt


def plot_result(iter_num):
    data_1 = np.load(f'Figure/Sol_iter{iter_num-1}_result.npz')
    data_2 = np.load(f'Figure/Err_iter{iter_num}_result.npz')

    t = data_1['t']
    h = data_1['h']
    v = data_1['v']
    m = data_1['mass']
    u = data_1['u']
    n_fe_list = data_1['n_fe_list']
    n_cp_list = data_1['n_cp_list']
    n_itv = len(n_fe_list)
    col_list = data_2['col_list']
    ncol_list = data_2['ncol_list']
    d_list = data_2['d_list']
    h_list = data_2['h_list']
    c_list = data_2['c_list']
    u_0 = data_2['u_0']

    plt.figure(1)
    plt.plot(t, h, color='#9400D3')
    plt.plot(t, h, '#9400D3', marker='o', markersize=4)
    plt.xlabel('t(s)')
    plt.ylabel('h')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_1.png')

    plt.figure(2)
    plt.plot(t, v, color='#9400D3')
    plt.plot(t, v, '#9400D3', marker='o', markersize=4)
    plt.xlabel('t(s)')
    plt.ylabel('v(m/s)')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_2.png')

    plt.figure(3)
    plt.plot(t, m, color='#9400D3')
    plt.plot(t, m, '#9400D3', marker='o', markersize=4)
    plt.xlabel('t(s)')
    plt.ylabel('m(kg)')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_3.png')

    plt.figure(4)
    u[0] = u[1]
    # Plot interval and finite-element
    count = 0
    plt.plot(t[0], u[0], '#9400D3', marker='s', markersize=4)
    for itv in range(n_itv):
        n_fe = n_fe_list[itv]
        n_cp = n_cp_list[itv]
        for fe in range(n_fe):
            for cp in range(n_cp):
                p = count + fe * n_cp + cp + 1
                if cp == n_cp - 1:
                    plt.plot(t[p], u[p], '#9400D3', marker='s', markersize=4)
                else:
                    plt.plot(t[p], u[p], '#9400D3', marker='^', markersize=4)
        count = count + n_fe * n_cp
        if itv != n_itv - 1:
            plt.plot(t[count], u[count], 'b', marker='s', markersize=4)
    plt.plot(t, u, color='#9400D3')
    plt.xlabel('t(s)')
    plt.ylabel('u')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_4.png')

    plt.figure(5)
    plt.plot(ncol_list, d_list, color='#9400D3')
    plt.plot(ncol_list, d_list, '#9400D3', marker='o', markersize=4)
    plt.xlabel('t(s)')
    plt.ylabel('defection')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_5.png')

    plt.figure(6)
    plt.plot(col_list, h_list, color='#9400D3')
    plt.plot(col_list, h_list, '#9400D3', marker='o', markersize=4)
    plt.xlabel('t(s)')
    plt.ylabel('hamilton')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_6.png')

    plt.figure(7)
    plt.plot(ncol_list, c_list, color='#9400D3')
    plt.plot(ncol_list, c_list, '#9400D3', marker='o', markersize=4)
    plt.xlabel('t(s)')
    plt.ylabel('constraint')
    plt.grid(ls='--')
    plt.tight_layout()
    plt.savefig(f'Figure/Figure_7.png')
    plt.show()
