from tardis.continuum.base import PhysicalContinuumProcess
from tardis.continuum.radiative_processes import *
from tardis.continuum.collisional_processes import *
from tardis.continuum.input_data import ContinuumInputData
from tardis.continuum.base import InverseProcess
from tardis.continuum.cooling import CoolingRates
from tardis.continuum.exceptions import InvalidContinuumProcessError
from tardis.continuum.probabilities import RecombinationTransitionProbabilities, TransitionProbabilities

default_processes = ['radiative_recombination', 'free_free', 'collisional_recombination', 'radiative_excitation',
                     'collisional_ionization', 'collisional_excitation', 'collisional_deexcitation',
                     'radiative_deexcitation', 'radiative_ionization']


class BaseContinuum(object):
    direct_processes = {process.name: process for process in PhysicalContinuumProcess.__subclasses__()}
    inverse_processes = [process.name for process in InverseProcess.__subclasses__()]
    process2inverse_process = {process.name_of_inverse_process: process for process in InverseProcess.__subclasses__()}

    def __init__(self, atom_data, plasma_array, ws, radiative_transition_probabilities, estimators,
                 requested_processes=default_processes):
        self._validate_requested_processes(requested_processes)
        self.input = ContinuumInputData(atom_data, plasma_array, ws, radiative_transition_probabilities, estimators)
        self._set_physical_processes(requested_processes)
        self._set_inverse_processes()
        self._set_cooling_rates()
        self._set_recombination_transition_probabilities()
        self._set_transition_probabilities()

    def _set_physical_processes(self, requested_processes):
        for name, process in BaseContinuum.direct_processes.iteritems():
            if name in requested_processes:
                setattr(self, name, process(self.input))

    def _set_inverse_processes(self):
        for name_of_inverse_process, process in BaseContinuum.process2inverse_process.iteritems():
            if hasattr(self, name_of_inverse_process):
                inverse_process = getattr(self, name_of_inverse_process)
                setattr(self, process.name, process.from_inverse_process(inverse_process))

    def _set_cooling_rates(self):
        cooling_processes = {name: process for name, process in self.__dict__.iteritems()
                             if hasattr(process, 'cooling') and process.cooling is True}
        self.cooling_rates = CoolingRates(self.input, **cooling_processes)

    def _set_transition_probabilities(self):
        macro_atom_processes = {name: process for name, process in self.__dict__.iteritems()
                                if hasattr(process, 'macro_atom_transitions') and
                                process.macro_atom_transitions is not None}
        self.transition_probabilities = TransitionProbabilities(self.input, **macro_atom_processes)


    def _set_recombination_transition_probabilities(self):
        recombination_process_list = ['radiative_recombination', 'collisional_recombination']
        recombination_processes = {}
        for process in recombination_process_list:
            if hasattr(self, process):
                recombination_processes[process] = getattr(self, process)
        self.recombination_transition_probabilities = \
            RecombinationTransitionProbabilities(self.input, **recombination_processes)

    @classmethod
    def _validate_requested_processes(cls, requested_processes):
        for process_name in requested_processes:
            if process_name in cls.direct_processes or process_name in cls.inverse_processes:
                continue
            else:
                raise InvalidContinuumProcessError(process_name)