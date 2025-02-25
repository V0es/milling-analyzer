from .mill import Mill


class Space:
    def __init__(
            self,
            mill_cutter: Mill,
            ):
        self.mill = mill_cutter


