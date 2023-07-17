import numpy as np
from tardis.model.geometry.cartesian3d import Cartesian3DGeometry


def cubegeometry(velocity, numbercells, timeofexplosion):
    x = np.linspace(-velocity, velocity, numbercells + 1)
    x_inner = x[:-1]
    x_outer = x[1:]
    v_inner = np.dstack(np.meshgrid([x_inner] * 3)).reshape(-1, 3)
    v_outer = np.dstack(np.meshgrid([x_outer] * 3)).reshape(-1, 3)
    # v = np.dstack(np.meshgrid([x]*3)).reshape(-1, 3)
    # r = v * timeofexplosion
    r_inner = v_inner * timeofexplosion
    r_outer = v_outer * timeofexplosion
    return Cartesian3DGeometry(
        r_inner, r_outer, v_inner, v_outer, [numbercells] * 3
    )
