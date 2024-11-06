from enum import Enum
from typing import List

import numpy as np

class MillingDirection(Enum):
    UP = 0
    DOWN = 1


class CutterTooth:
    def __init__(self, initial_angle: float):
        self.initial_angle = initial_angle

    def angle(self, rpm: float, time: float) -> float:
        return (2 * np.pi * rpm / 60) * time + self.initial_angle


class Mill:
    def __init__(
            self,
            teeth_num: int,
            tangential_force_coeff: float,
            normal_force_coeff: float,
            relative_damping: float,
            aD: float,
            natural_frequency: float = None,
            modal_mass: float = None,
            direction: MillingDirection = MillingDirection.DOWN
            ):
        self.teeth_num = teeth_num
        self.tangential_force_coeff = tangential_force_coeff
        self.normal_force_coeff = normal_force_coeff
        self.angular_natural_frequency = natural_frequency * 2 * np.pi
        self.relative_damping = relative_damping
        self.modal_mass = modal_mass
        self.aD = aD
        self.direction = direction
        self.teeth: List[CutterTooth] = []
        if direction == MillingDirection.DOWN:
            self.start_angle = np.arccos(2 * aD - 1)
            self.exit_angle = np.pi
        else:
            self.start_angle = 0
            self.exit_angle = np.arccos(2 * aD - 1)

        self.initialize_tooth()

    def initialize_tooth(self):
        for i in range(self.teeth_num):
            angle = 2 * np.pi * i / self.teeth_num
            tooth = CutterTooth(angle)
            self.teeth.append(tooth)

    def tooth_in_cut(self, tooth_angle: float) -> bool:
        if self.start_angle < tooth_angle < self.exit_angle:
            return True
        else:
            return False

