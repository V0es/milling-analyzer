from dataclasses import dataclass
from enum import Enum, StrEnum
from typing import Callable, Tuple, List

import numpy as np
from numpy.linalg import inv, eigvals
from scipy.linalg import expm

from src.mill import Mill

class VariableName(StrEnum):
    """
    Enum class for variable parameters
    """
    SPINDLE_SPEED = 'Spindel Speeds, m'
    ANGULAR_NATURAL_FREQUENCY = 'Natural Frequencies, rad/s'
    DEPTH_OF_CUT = 'Depths Of Cut, m'
    MODAL_MASS = 'Modal Masses, kg'


@dataclass
class Varying:
    """
    Dataclass to hold information about name of variable and calculating parameters
    """
    var_name: str
    start_value: float
    final_value: float
    steps: int


class Solver:
    def __init__(
            self,
            mill: Mill,
            x_variable: Varying,
            y_variable: Varying,
            cutter_function: Callable = None,
            intervals_per_period: int = 40,
            integration_steps: int = 20,
            weight_a: float = 0.5,
            weight_b: float = 0.5,


    ):
        self.x_var = x_variable
        self.y_var = y_variable
        self.mill_cutter = mill

        self.integration_steps = integration_steps
        self.weight_a = weight_a
        self.weight_b = weight_b
        if cutter_function:
            self.vec_cutter_func = np.vectorize(cutter_function)
        self.intervals_per_period = intervals_per_period

        self.D = np.zeros((self.intervals_per_period+2, self.intervals_per_period+2))
        self.d = np.ones((self.intervals_per_period+1, 1))
        self.d[0:3] = 0
        self.D += np.diag(self.d, -1)
        self.D[2, 0] = 1

    def integrate_force_function(self) -> np.ndarray:
        h_i = np.zeros(self.intervals_per_period)
        for i in range(self.intervals_per_period):
            dtr = 2 * np.pi / (self.mill_cutter.teeth_num * self.intervals_per_period)
            for j in range(self.mill_cutter.teeth_num):
                for h in range(self.integration_steps):
                    fi = i * dtr + j * 2 * np.pi / self.mill_cutter.teeth_num + h * dtr / self.integration_steps
                    h_i[i] +=\
                        (self.mill_cutter.tooth_in_cut(fi) *
                         (self.mill_cutter.tangential_force_coeff * np.cos(fi) +
                          self.mill_cutter.normal_force_coeff * np.sin(fi)) * np.sin(fi) / self.integration_steps)
        return h_i

    def solve(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

        ss = np.zeros((self.x_var.steps, self.y_var.steps))
        dc = np.zeros((self.x_var.steps, self.y_var.steps))
        ei = np.zeros((self.x_var.steps, self.y_var.steps))

        h_i = self.integrate_force_function()

        D = np.zeros([self.intervals_per_period + 2, self.intervals_per_period + 2])
        d = np.ones([self.intervals_per_period + 1])
        d[0:2] = 0
        D += np.diag(d, -1)
        D[2][0] = 1

        for x in range(self.x_var.steps):
            o = self.x_var.start_value + x * (self.x_var.final_value - self.x_var.start_value) / self.x_var.steps
            tau = 60 / (o * self.mill_cutter.teeth_num)
            dt = tau / self.intervals_per_period

            for y in range(self.y_var.steps):
                w = self.y_var.start_value + y * (self.y_var.final_value - self.y_var.start_value) / self.y_var.steps
                Fi = np.eye(self.intervals_per_period + 2)
                for i in range(self.intervals_per_period):
                    A = np.zeros([2, 2])
                    A[0, 1] = 1
                    A[1, 0] = -self.mill_cutter.angular_natural_frequency ** 2 - h_i[i] * w / self.mill_cutter.modal_mass
                    A[1, 1] = -2 * self.mill_cutter.relative_damping * self.mill_cutter.angular_natural_frequency
                    B = np.zeros((2, 2))
                    B[1, 0] = h_i[i] * w / self.mill_cutter.modal_mass
                    P = expm(A * dt)

                    R = (expm(A * dt) - np.eye(2)).dot(inv(A).dot(B))
                    D[:2, :2] = P
                    D[:2, self.intervals_per_period] = self.weight_a * R[:, 0]
                    D[0:2, self.intervals_per_period + 1] = self.weight_b * R[:, 0]
                    Fi = np.dot(D, Fi)
                ss[x, y] = o
                dc[x, y] = w
                ei[x, y] = max(abs(eigvals(Fi)))
            print(self.x_var.steps + 1 - x)

        return ss, dc, ei

