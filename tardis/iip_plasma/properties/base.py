import logging
from abc import ABCMeta, abstractmethod, abstractproperty

import numpy as np
import pandas as pd

__all__ = [
    "ArrayInput",
    "BaseAtomicDataProperty",
    "BasePlasmaProperty",
    "ConvergedPlasmaProperty",
    "DataFrameInput",
    "HiddenPlasmaProperty",
    "Input",
    "PreviousIterationProperty",
    "ProcessingPlasmaProperty",
]

logger = logging.getLogger(__name__)


class BasePlasmaProperty:
    """
    Attributes
    ----------
    outputs : Tuple (strings)
              List of output parameter names for particular property.
    name : String
           Class name
    latex_name : String
                 Used to label nodes when plotting graphs
    """

    __metaclass__ = ABCMeta

    @abstractproperty
    def outputs(self):
        pass

    @property
    def name(self):
        return self.__class__.__name__

    def __init__(self):
        for output in self.outputs:
            setattr(self, output, None)

    def _update_type_str(self):
        self.type_str = repr(type(self.value))

    def get_latex_label(self):
        latex_template = r"""\textbf{{Name}} {name}
\textbf{{Formula}} {formula}
{description}
"""
        outputs = self.outputs.replace("_", r"\_")
        latex_name = getattr(self, "latex_name", "")
        if latex_name != "":
            complete_name = f"{latex_name} [{self.latex_name}]"
        else:
            complete_name = latex_name

        latex_label = latex_template.format(
            name=complete_name,
            formula=getattr(self, "latex_formula", "--"),
            description=getattr(self, "latex_description", ""),
        )
        return latex_label.replace("\\", r"\\")


class ProcessingPlasmaProperty(BasePlasmaProperty):
    """
    Attributes
    ----------
    inputs : Tuple (strings)
             List of input parameters required to create the property
    """

    __metaclass__ = ABCMeta

    def __init__(self, plasma_parent):
        super(ProcessingPlasmaProperty, self).__init__()
        self.plasma_parent = plasma_parent
        self._update_inputs()

    def _update_inputs(self):
        """
        This function uses the CPython API to read the variable names from the
        `calculate`-function and makes the plasma routines easily programmable.
        """
        calculate_call_signature = self.calculate.__code__.co_varnames[
            : self.calculate.__code__.co_argcount
        ]
        self.inputs = [
            item for item in calculate_call_signature if item != "self"
        ]

    def _get_input_values(self):
        return (self.plasma_parent.get_value(item) for item in self.inputs)

    def update(self):
        """
        Updates the processing Plasma by calling the `calculate`-method with
        the required inputs

        :return:
        """
        if len(self.outputs) == 1:
            setattr(
                self, self.outputs[0], self.calculate(*self._get_input_values())
            )
        else:
            new_values = self.calculate(*self._get_input_values())
            for i, output in enumerate(self.outputs):
                setattr(self, output, new_values[i])

    @abstractmethod
    def calculate(self, *args, **kwargs):
        raise NotImplementedError(
            "This method needs to be implemented by processing plasma modules"
        )


class HiddenPlasmaProperty(ProcessingPlasmaProperty):
    """
    Used for plasma properties that should not be displayed in the final graph (e.g. lines_lower_level_index).
    The code will automatically remove these property names from the graph and instead connect their inputs directly
    to their outputs.
    """

    __metaclass__ = ABCMeta

    def __init__(self, plasma_parent):
        super(HiddenPlasmaProperty, self).__init__(plasma_parent)


class BaseAtomicDataProperty(ProcessingPlasmaProperty):
    """
    Used for atomic data properties. Main feature is the ability to filter atomic data by the elements required for
    the simulation.
    """

    __metaclass__ = ABCMeta

    inputs = ["atomic_data", "selected_atoms"]

    def __init__(self, plasma_parent):
        super(BaseAtomicDataProperty, self).__init__(plasma_parent)
        self.value = None

    @abstractmethod
    def _set_index(self, raw_atomic_property):
        raise NotImplementedError("Needs to be implemented in subclasses")

    @abstractmethod
    def _filter_atomic_property(self, raw_atomic_property):
        raise NotImplementedError("Needs to be implemented in subclasses")

    def calculate(self, atomic_data, selected_atoms):
        if getattr(self, self.outputs[0]) is not None:
            return getattr(self, self.outputs[0])
        # Atomic Data Issue: Some atomic property names in the h5 files are preceded
        # by an underscore, e.g. _levels, _lines.
        try:
            raw_atomic_property = getattr(atomic_data, "_" + self.outputs[0])
        except AttributeError:
            raw_atomic_property = getattr(atomic_data, self.outputs[0])
        finally:
            return self._set_index(
                self._filter_atomic_property(
                    raw_atomic_property, selected_atoms
                )
            )


class Input(BasePlasmaProperty):
    """
    The plasma property class for properties that are input directly from model and not calculated within the plasma
    module, e.g. t_rad.
    """

    def _set_output_value(self, output, value):
        setattr(self, output, value)

    def set_value(self, value):
        assert len(self.outputs) == 1
        self._set_output_value(self.outputs[0], value)


class ArrayInput(Input):
    def _set_output_value(self, output, value):
        setattr(self, output, np.array(value, copy=False))


class DataFrameInput(Input):
    def _set_output_value(self, output, value):
        setattr(self, output, np.array(pd.DataFrame(value), copy=False))


class PreviousIterationProperty(BasePlasmaProperty):
    """
    This class is used for properties where, to prevent a property calculation loop, the property values from the
    previous iteration (which are static) are used in the current calculation. Usually only required for NLTE
    calculations. Given a sufficient number of iterations, the values should converge successfully on the correct
    solution.
    """

    def _set_initial_value(self, value):
        self.set_value(value)

    def _set_output_value(self, output, value):
        setattr(self, output, value)

    def set_value(self, value):
        assert len(self.outputs) == 1
        self._set_output_value(self.outputs[0], value)


class ConvergedPlasmaProperty(ProcessingPlasmaProperty):
    """
    This class is used for properties which do not need to be calculated
    during an iterative NLTE calculation.
    """
