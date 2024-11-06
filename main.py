import time

from matplotlib import pyplot as plt

from src.mill import Mill, MillingDirection
from src.solver import Solver, Varying, VariableName


mill = Mill(
    teeth_num=2,
    tangential_force_coeff=6e8,
    normal_force_coeff=2e8,
    relative_damping=0.011,
    aD=0.05,
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

if __name__ == '__main__':
    t0 = time.time()
    ss, dc, ei = solver.solve()
    print(f"Время выполнения: {time.time() - t0} секунд")
    plt.figure()
    plt.contour(ss, dc, ei, [1], colors='k')
    plt.xlabel('Spindle Speed (rpm)')
    plt.ylabel('Depth of Cut (m)')
    plt.title('Stability Contour')
    plt.show()