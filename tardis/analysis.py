#codes to for analyse the model.

from astropy import units as u
import numpy as np


def get_last_line_interaction(wavelength_start, wavelength_end, model):
    nu_start = wavelength_end.to('Hz', u.spectral())
    nu_end = wavelength_start.to('Hz', u.spectral())
    wavelength_filter = (model.montecarlo_nu > nu_start) & (model.montecarlo_nu < nu_end) & \
                        (model.last_line_interaction_in_id != -1)

    last_line_in_ids = model.last_line_interaction_in_id[wavelength_filter]
    last_line_out_ids = model.last_line_interaction_out_id[wavelength_filter]

    line_in_table = model.atom_data.lines.ix[last_line_in_ids].groupby(['atomic_number', 'ion_number'])[
        'wavelength'].count()
    line_in_table = line_in_table.astype(float) / line_in_table.sum()

    line_out_table = model.atom_data.lines.ix[last_line_out_ids].groupby(['atomic_number', 'ion_number'])[
        'wavelength'].count()
    line_out_table = line_out_table.astype(float) / line_out_table.sum()

    print "Lines in:\n------\n %s" % line_in_table
    print "Lines out:\n------\n %s" % line_out_table

    return last_line_in_ids, last_line_out_ids


class LastLineInteraction(object):
    def __init__(self, model):
        self.model = model
        self._wavelength_start = 0 * u.m
        self._wavelength_end = 1 * u.m
        self._atomic_number = None
        self._ion_number = None
        self.update_last_interaction()

    @property
    def wavelength_start(self):
        return self._wavelength_start.to('angstrom')

    @wavelength_start.setter
    def wavelength_start(self, value):
        if not isinstance(value, u.Quantity):
            raise ValueError('needs to be a Quantity')
        self._wavelength_start = value
        self.update_last_interaction()

    @property
    def wavelength_end(self):
        return self._wavelength_end.to('angstrom')

    @wavelength_end.setter
    def wavelength_end(self, value):
        if not isinstance(value, u.Quantity):
            raise ValueError('needs to be a Quantity')
        self._wavelength_end = value
        self.update_last_interaction()

    @property
    def atomic_number(self):
        return self._atomic_number

    @atomic_number.setter
    def atomic_number(self, value):
        self._atomic_number = value
        self.update_last_interaction()

    @property
    def ion_number(self):
        return self._ion_number

    @ion_number.setter
    def ion_number(self, value):
        self._ion_number = value
        self.update_last_interaction()

    def update_last_interaction(self):
        last_line_interaction_in_id = self.model.last_line_interaction_in_id[
            self.model.last_line_interaction_in_id != -1]
        last_line_interaction_out_id = self.model.last_line_interaction_out_id[
            self.model.last_line_interaction_out_id != -1]

        line_list_in = self.model.atom_data.lines.ix[last_line_interaction_in_id]
        line_list_out = self.model.atom_data.lines.ix[last_line_interaction_out_id]

        mask_in = (line_list_out['wavelength'] > self.wavelength_start) & (
        line_list_out['wavelength'] < self.wavelength_end)
        mask_out = mask_in.copy()
        if self.atomic_number is not None:
            mask_in &= line_list_in.atomic_number == self.atomic_number
            mask_out &= line_list_out.atomic_number == self.atomic_number

        if self.ion_number is not None:
            mask_in &= line_list_in.ion_number == self.ion_number
            mask_out &= line_list_out.ion_number == self.ion_number

        self.last_line_list_in = line_list_in[mask_in]
        self.last_line_list_out = line_list_out[mask_out]

        self.contributing_last = self.last_line_list_in.groupby(['atomic_number', 'ion_number'])['wavelength'].count()
        self.contributing_last = self.contributing_last.astype(float) / self.contributing_last.sum()





