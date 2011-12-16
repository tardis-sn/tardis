#File contains model class
class Model(object):
    pass

class OneZoneBaseModel(Model):
    def __init__(self, r1, r2, v1):
        """Initialize the model
        Parameters
        r1 : float
            inner radius
        r2 : float
            outer radius
        v1 : float
            velocity at inner radius