from dataclasses import dataclass

import numpy as np

from tardis.plasma.radiation_field import DilutePlanckianRadiationField


@dataclass
class EstimatedRadiationFieldProperties:
    dilute_blackbody_radiationfield_state: DilutePlanckianRadiationField
    j_blues: np.ndarray
