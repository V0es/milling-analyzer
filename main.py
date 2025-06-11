import time

import numpy as np
from matplotlib import pyplot as plt

from src.mill import Mill, MillingDirection
from src.solver import Solver, Varying, VariableName
from src.solver2dof import Solver2DOF

X_STEPS = 100
Y_STEPS = 50
BASE_FREQ = 922
BASE_DAMP = 0.011
BASE_AD = 0.05
TEETH_NUM = 2
MODAL_MASS = 0.03993

omega_st_factor = 5.532
omega_fin_factor = 1.1064


def get_1dof_damp_misl(damps):
    ei_new_list = []
    misls = []
    ei_base = None
    for damp in damps:
        x_var = Varying(
            var_name=VariableName.SPINDLE_SPEED,
            start_value=5e3,
            final_value=25e3,
            steps=X_STEPS
        )

        y_var = Varying(
            var_name=VariableName.DEPTH_OF_CUT,
            start_value=0e-3,
            final_value=35e-3,
            steps=Y_STEPS
        )
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=damp,
            aD=BASE_AD,
            natural_frequency=BASE_FREQ,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )
        solver = Solver(
            mill,
            x_variable=x_var,
            y_variable=y_var
        )
        ss, dc, ei = solver.solve_jit_test()
        if damp == BASE_DAMP:
            ei_base = ei
        ei_new_list.append(ei)
    for ei in ei_new_list:
        misls.append(get_misl(ei_base, ei, dc))
    print('damp 1dof misls', misls)
    return misls

def get_1dof_freq_misl(freqs):
    ei_new_list = []
    misls = []
    ei_base = None
    for freq in freqs:
        om_st = 60 * freq / (2 * omega_st_factor)
        om_fin = 60 * freq / (2 * omega_fin_factor)
        x_var = Varying(
            var_name=VariableName.SPINDLE_SPEED,
            start_value=om_st,
            final_value=om_fin,
            steps=X_STEPS
        )

        y_var = Varying(
            var_name=VariableName.DEPTH_OF_CUT,
            start_value=0e-3,
            final_value=35e-3,
            steps=Y_STEPS
        )
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=BASE_DAMP,
            aD=BASE_AD,
            natural_frequency=freq,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )
        solver = Solver(
            mill,
            x_variable=x_var,
            y_variable=y_var
        )
        ss, dc, ei = solver.solve_jit_test()
        if freq == BASE_FREQ:
            ei_base = ei
        ei_new_list.append(ei)
    for ei in ei_new_list:
        misls.append(get_misl(ei_base, ei, dc))

    return misls


def plot_aD_compare_misl(ads, ei_base_list, ei_new_list, depths):
    fig, ax = plt.subplots()

    categories = [f'a/D = {ad}' for ad in ads]
    misls = []

    for ei_base, ei_new in zip(ei_base_list, ei_new_list):
        misls.append(get_misl(ei_base, ei_new, depths))

    ax.bar(categories, misls)
    ax.grid()
    ax.set_xlabel('Радиальная глубина резания a/D')
    ax.set_ylabel('Среднее увеличение предела устойчивости, м')
    fig.suptitle('Среднее увеличение предела устойчивости для разных значений a/D относительно 1-DOF модели')
    fig.tight_layout()


def plot_freq_compare_misl(freqs, ei_base, ei_new_list, depths):
    fig, ax = plt.subplots()
    categories = [fr'$f_n$ = {round(freq, 1)} Гц' for freq in freqs]
    misls = []

    for ei_new in ei_new_list:
        misls.append(get_misl(ei_base, ei_new, depths))

    ax.bar(categories, misls)
    ax.grid()
    ax.set_xlabel(r'Cобственная частота $f_n$, Гц')
    ax.set_ylabel('Среднее увеличение предела устойчивости, м')
    fig.suptitle(r'Среднее увеличение предела устойчивости для разных значений $f_n$')
    fig.tight_layout()

    plot_grouped_bars_freq_misl(freqs, misls)

def plot_damp_compare_misl(damps, ei_base, ei_new_list, depths):
    fig, ax = plt.subplots()
    categories = [fr'$\xi$ = {round(damp, 5)}' for damp in damps]
    misls = []

    for ei_new in ei_new_list:
        misls.append(get_misl(ei_base, ei_new, depths))

    ax.bar(categories, misls)
    ax.grid()
    ax.set_xlabel(r'Относительное демпфирование $\xi$')
    ax.set_ylabel('Среднее увеличение предела устойчивости, м')
    fig.suptitle(r'Среднее увеличение предела устойчивости для разных значений $\xi$')
    fig.tight_layout()

    plot_grouped_bars_damp_misl(damps, misls)


def plot_damps_compare(damps, x_var, y_var):

    fig, ax = plt.subplots()
    alphas = np.linspace(1, 0.2, len(damps))
    # axes = np.ravel(axes)
    legend_elements = []
    legend_labels = []
    ei_new_list = []
    ei_base = None
    for damp, alpha in zip(damps, alphas):

        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=damp,
            aD=BASE_AD,
            natural_frequency=BASE_FREQ,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )
        solver = Solver2DOF(
            mill,
            x_variable=x_var,
            y_variable=y_var
        )
        ss, dc, ei = solver.solve_jit_test()
        if damp == BASE_DAMP:
            ei_base = ei
        ei_new_list.append(ei)
        #TODO: убрать костыль
        if damp < BASE_DAMP*0.7 or damp > BASE_DAMP*1.3:
            continue

        cp_freq = ax.contour((mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss * mill.teeth_num), dc, ei, [1], colors='k', alpha=alpha)
        # ax.set_title(fr'$f_n$ = {round(freq, 1)} Гц')
        leg_el, _ = cp_freq.legend_elements()
        legend_elements.append(leg_el[0])
        legend_labels.append(rf'$\xi$ = {round(damp, 5)}')
        ax.set_xlabel(r'$\omega_f$')
        ax.set_ylabel('Глубина резания, м')
    fig.suptitle(fr'Диаграммы устойчивости для разных значений $\xi$, $f_n$ = {BASE_FREQ} Гц, a/D = {BASE_AD}')
    ax.legend(legend_elements, legend_labels)
    fig.tight_layout()

    # plot_damp_compare_misl(damps, ei_base, ei_new_list, dc)


def plot_freqs_compare(freqs):
    fig, ax = plt.subplots()
    alphas = np.linspace(1, 0.2, len(freqs))
    # axes = np.ravel(axes)
    legend_elements = []
    legend_labels = [rf'$f_n$ = {round(freq, 1)} Гц' for freq in freqs]
    ei_new_list = []
    ei_base = None
    for freq, alpha in zip(freqs, alphas):
        om_st = 60 * freq / (2 * omega_st_factor)
        om_fin = 60 * freq / (2 * omega_fin_factor)
        x_var = Varying(
            var_name=VariableName.SPINDLE_SPEED,
            start_value=om_st,
            final_value=om_fin,
            steps=X_STEPS
        )

        y_var = Varying(
            var_name=VariableName.DEPTH_OF_CUT,
            start_value=0e-3,
            final_value=35e-3,
            steps=Y_STEPS
        )
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=BASE_DAMP,
            aD=BASE_AD,
            natural_frequency=freq,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )
        solver = Solver2DOF(
            mill,
            x_variable=x_var,
            y_variable=y_var
        )
        ss, dc, ei = solver.solve_jit_test()
        if freq == BASE_FREQ:
            ei_base = ei
        ei_new_list.append(ei)
        cp_freq = ax.contour((mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss * mill.teeth_num), dc, ei, [1], colors='k', alpha=alpha)
        # ax.set_title(fr'$f_n$ = {round(freq, 1)} Гц')
        leg_el, _ = cp_freq.legend_elements()
        legend_elements.append(leg_el[0])
        ax.set_xlabel(r'$\omega_f$')
        ax.set_ylabel('Глубина резания, м')
    fig.suptitle(fr'Диаграммы устойчивости для разных значений $f_n$, $\xi$ = {BASE_DAMP}, a/D = {BASE_AD}')
    ax.legend(legend_elements, legend_labels)
    fig.tight_layout()

    plot_freq_compare_misl(freqs, ei_base, ei_new_list, dc)


def plot_freqs_separate(freqs):
    fig, axes = plt.subplots(2, 2)
    axes = np.ravel(axes)
    for ax, freq in zip(axes, freqs):
        om_st = 60 * freq / (2 * omega_st_factor)
        om_fin = 60 * freq / (2 * omega_fin_factor)
        x_var = Varying(
            var_name=VariableName.SPINDLE_SPEED,
            start_value=om_st,
            final_value=om_fin,
            steps=X_STEPS
        )

        y_var = Varying(
            var_name=VariableName.DEPTH_OF_CUT,
            start_value=0e-3,
            final_value=35e-3,
            steps=Y_STEPS
        )
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=BASE_DAMP,
            aD=BASE_AD,
            natural_frequency=freq,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )
        solver = Solver2DOF(
            mill,
            x_variable=x_var,
            y_variable=y_var
        )

        ss, dc, ei = solver.solve_jit_test()
        ax.contour((mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss * mill.teeth_num), dc, ei, [1], colors='k')
        ax.set_title(fr'$f_n$ = {round(freq, 1)} Гц')
        ax.set_xlabel(r'$\omega_f$')
        ax.set_ylabel('Глубина резания, м')
    fig.suptitle(fr'Диаграммы устойчивости для разных значений $f_n$, $\xi$ = {BASE_DAMP}, a/D = {BASE_AD}')
    fig.tight_layout()


def plot_damps_separate(damps, x_var, y_var):
    fig, axes = plt.subplots(2, 2)
    axes = np.ravel(axes)
    for ax, damp in zip(axes, damps):
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=damp,
            aD=BASE_AD,
            natural_frequency=BASE_FREQ,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )
        solver = Solver2DOF(
            mill,
            x_variable=x_var,
            y_variable=y_var
        )

        ss, dc, ei = solver.solve_jit_test()
        ax.contour((mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss * mill.teeth_num), dc, ei, [1], colors='k')
        ax.set_title(fr'$\xi$ = {round(damp, 5)}')
        ax.set_xlabel(r'$\omega_f$')
        ax.set_ylabel('Глубина резания, м')
    fig.suptitle(fr'Диаграммы устойчивости для разных значений $\xi$, $f_n$ = {BASE_FREQ} Гц, a/D = {BASE_AD}')
    fig.tight_layout()


def plot_ads_compare(ads, x_var, y_var):
    fig, axes = plt.subplots(4, 4)
    axes = np.ravel(axes)

    ei_1dof_list = []
    ei_2dof_list = []
    for ad, ax in zip(ads, axes):
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=BASE_DAMP,
            aD=ad,
            natural_frequency=BASE_FREQ,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )

        solver1d = Solver(
            mill=mill,
            x_variable=x_var,
            y_variable=y_var
        )

        solver2d = Solver2DOF(
            mill=mill,
            x_variable=x_var,
            y_variable=y_var
        )

        ss_1d, dc_1d, ei_1d = solver1d.solve_jit_test()
        ss_2d, dc_2d, ei_2d = solver2d.solve_jit_test()

        ei_1dof_list.append(ei_1d)
        ei_2dof_list.append(ei_2d)

        cp_1dof = ax.contour(
            (mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss_1d * mill.teeth_num),
            dc_1d,
            ei_1d,
            [1],
            colors='k',
            label='1-DOF')

        cp_2dof = ax.contour(
            (mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss_2d * mill.teeth_num),
            dc_2d,
            ei_2d,
            [1],
            colors='k',
            linestyles='dashed',
            label='2-DOF')

        ax.set_title(f'a/D = {ad}')
        ax.set_xlabel(r'$\omega_f$')
        ax.set_ylabel('Глубина резания, м')

        h1, _ = cp_1dof.legend_elements()
        h2, _ = cp_2dof.legend_elements()
        ax.legend([h1[0], h2[0]], ['1-DOF', '2-DOF'])

    fig.suptitle(fr'Диаграммы устойчивости для разных значений a/D, $\xi$ = {BASE_DAMP}, $f_n$ = {BASE_FREQ} Гц')
    fig.tight_layout()

    plot_aD_compare_misl(ads, ei_1dof_list, ei_2dof_list, dc_1d)

def plot_grouped_bars_damp_misl(damps, misl_2dof):
    bar_width = 0.3
    misl_1dof = get_1dof_damp_misl(damps)
    categories = [fr'$\xi$ = {round(damp, 5)}' for damp in damps]
    fig, ax = plt.subplots()
    x = np.arange(len(categories))
    ax.bar(x, misl_1dof, width=bar_width, label='1-DOF', color='gray')
    ax.bar(x + bar_width, misl_2dof, width=bar_width, label='2-DOF', color='black')

    ax.set_xticks(x, categories)
    ax.grid()
    ax.set_xlabel(r'Относительное демпфирование $\xi$')
    ax.set_ylabel('Среднее увеличение предела устойчивости, м')
    fig.suptitle(r'Среднее увеличение предела устойчивости для разных значений $\xi$: сравнение 1-DOF и 2-DOF')
    ax.legend()
    fig.tight_layout()



def plot_grouped_bars_freq_misl(freqs, misl_2dof):
    bar_width = 0.3
    misl_1dof = get_1dof_freq_misl(freqs)
    categories = [fr'$f_n$ = {round(freq, 1)} Гц' for freq in freqs]
    fig, ax = plt.subplots()
    x = np.arange(len(categories))
    ax.bar(x, misl_1dof, width=bar_width, label='1-DOF', color='gray')
    ax.bar(x + bar_width, misl_2dof, width=bar_width, label='2-DOF', color='black')

    ax.set_xticks(x, categories)
    ax.grid()
    ax.set_xlabel(r'Cобственная частота $f_n$, Гц')
    ax.set_ylabel('Среднее увеличение предела устойчивости, м')
    fig.suptitle(r'Среднее увеличение предела устойчивости для разных значений $f_n$: сравнение 1-DOF и 2-DOF')
    ax.legend()
    fig.tight_layout()


def plot_ads(ads, x_var, y_var):
    fig, axes = plt.subplots(2, 2)
    axes = np.ravel(axes)
    for ad, ax in zip(ads, axes):
        mill = Mill(
            teeth_num=TEETH_NUM,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=BASE_DAMP,
            aD=ad,
            natural_frequency=BASE_FREQ,
            modal_mass=MODAL_MASS,
            direction=MillingDirection.DOWN
        )

        solver = Solver2DOF(
            mill=mill,
            x_variable=x_var,
            y_variable=y_var
        )

        ss, dc, ei = solver.solve_jit_test()

        ax.contour((mill.angular_natural_frequency * 60 / (2 * np.pi)) / (ss * mill.teeth_num), dc, ei, [1], colors='k')
        ax.set_title(f'a/D = {ad}')
        ax.set_xlabel(r'$\omega_f$')
        ax.set_ylabel('Глубина резания, м')
    fig.suptitle(fr'Диаграммы устойчивости для разных значений a/D, $\xi$ = {BASE_DAMP}, $f_n$ = {BASE_FREQ} Гц')
    fig.legend()


def get_misl(ei_base: np.ndarray, ei_new: np.ndarray, depths: np.ndarray) -> float:
    ei_base_mask = np.where(ei_base <= 1, 1, 0)
    ei_new_mask = np.where(ei_new <= 1, 1, 0)

    base_depths = ei_base_mask * depths
    new_depths = ei_new_mask * depths

    base_critical_depths = np.max(base_depths, axis=0)
    new_critical_depths = np.max(new_depths, axis=0)

    misl = np.mean(new_critical_depths - base_critical_depths)

    return misl


def main():
    ads = np.linspace(0.01, 0.15, 9)

    freqs_sep = [BASE_FREQ * 0.7, BASE_FREQ * 0.85, BASE_FREQ * 1.15, BASE_FREQ * 1.3]
    freqs_comp = [BASE_FREQ*0.55, BASE_FREQ * 0.7, BASE_FREQ * 0.85, BASE_FREQ, BASE_FREQ * 1.15, BASE_FREQ * 1.3, BASE_FREQ*1.45]

    damps_sep = [BASE_DAMP * 0.7, BASE_DAMP * 0.85, BASE_DAMP * 1.15, BASE_DAMP * 1.3]
    damps_comp = [BASE_DAMP*0.55, BASE_DAMP * 0.7, BASE_DAMP * 0.85, BASE_DAMP, BASE_DAMP * 1.15, BASE_DAMP * 1.3, BASE_DAMP*1.45]

    x_var = Varying(
        var_name=VariableName.SPINDLE_SPEED,
        start_value=5e3,
        final_value=25e3,
        steps=X_STEPS
    )

    y_var = Varying(
        var_name=VariableName.DEPTH_OF_CUT,
        start_value=0e-3,
        final_value=35e-3,
        steps=Y_STEPS
    )
    plot_ads_compare(ads, x_var, y_var)
    # plot_ads(ads, x_var, y_var)
    # plot_freqs_separate(freqs_sep)
    # plt.tight_layout()
    # plot_freqs_compare(freqs)
    # plot_damps_separate(damps_sep, x_var, y_var)
    # plot_damps_compare(damps_comp, x_var, y_var)
    plt.show()


if __name__ == '__main__':
    main()












