import time

import numpy as np
from matplotlib import pyplot as plt

from src.mill import Mill, MillingDirection
from src.solver import Solver, Varying, VariableName

if __name__ == '__main__':
    ads = [0.1, 0.5, 1]
    base_damp = 0.011
    for idx, ad in enumerate(ads):
        mill = Mill(
            teeth_num=2,
            tangential_force_coeff=6e8,
            normal_force_coeff=2e8,
            relative_damping=base_damp,
            aD=ad,
            natural_frequency=922,
            modal_mass=0.03993,
            direction=MillingDirection.DOWN
        )

        x_var = Varying(
            var_name=VariableName.SPINDLE_SPEED,
            start_value=5e3,
            final_value=25e3,
            steps=400
        )

        y_var = Varying(
            var_name=VariableName.DEPTH_OF_CUT,
            start_value=0e-3,
            final_value=10e-3,
            steps=200
        )

        solver = Solver(
            mill=mill,
            x_variable=x_var,
            y_variable=y_var
        )

        t0 = time.time()
        ss, dc, ei = solver.solve_jit()
        print(f"Время выполнения: {time.time() - t0} секунд")
        plt.figure(idx)
        plot = plt.contour((mill.angular_natural_frequency / (2 * np.pi))*(ss*mill.teeth_num), dc, ei, [1])
        plt.clabel(plot, inline=1, fontsize=10)
        plt.xlabel(r'$\omega_f$')
        plt.ylabel('Depth of Cut (m)')
        plot.collections[0].set_label(fr'$\xi = {base_damp}$')
        plt.title(fr'Stability Contours, freq = {922}')
    plt.show()
