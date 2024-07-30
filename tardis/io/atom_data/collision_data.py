class CollisionData:
    def __init__(
        self,
        collision_energies=None,
        temperatures=None,
        thermally_averaged_strengths=None,
    ):
        self.data = collision_energies
        self.temperatures = temperatures
        self.yg = thermally_averaged_strengths

    def active_data(self):
        pass
