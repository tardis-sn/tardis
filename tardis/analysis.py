#codes to for analyse the model.

import re
import os

from astropy import units as u, constants
import numpy as np
import pandas as pd


class LastLineInteraction(object):

    @classmethod
    def from_model(cls, model):
        return cls(model.last_line_interaction_in_id, model.last_line_interaction_out_id,
                   model.last_line_interaction_shell_id, model.last_line_interaction_angstrom, model.atom_data.lines)

    def __init__(self, last_line_interaction_in_id, last_line_interaction_out_id, last_line_interaction_shell_id,
                 last_line_interaction_angstrom, lines, packet_filter_mode='packet_nu'):
        self.last_line_interaction_in_id = last_line_interaction_in_id
        self.last_line_interaction_out_id = last_line_interaction_out_id
        self.last_line_interaction_shell_id = last_line_interaction_shell_id
        self.last_line_interaction_angstrom = last_line_interaction_angstrom
        self.lines = lines

        self._wavelength_start = 0 * u.angstrom
        self._wavelength_end = np.inf * u.angstrom
        self._atomic_number = None
        self._ion_number = None
        self.packet_filter_mode = packet_filter_mode
        self.update_last_interaction_filter()



    @property
    def wavelength_start(self):
        return self._wavelength_start.to('angstrom')

    @wavelength_start.setter
    def wavelength_start(self, value):
        if not isinstance(value, u.Quantity):
            raise ValueError('needs to be a Quantity')
        self._wavelength_start = value
        self.update_last_interaction_filter()

    @property
    def wavelength_end(self):
        return self._wavelength_end.to('angstrom')

    @wavelength_end.setter
    def wavelength_end(self, value):
        if not isinstance(value, u.Quantity):
            raise ValueError('needs to be a Quantity')
        self._wavelength_end = value
        self.update_last_interaction_filter()

    @property
    def atomic_number(self):
        return self._atomic_number

    @atomic_number.setter
    def atomic_number(self, value):
        self._atomic_number = value
        self.update_last_interaction_filter()

    @property
    def ion_number(self):
        return self._ion_number

    @ion_number.setter
    def ion_number(self, value):
        self._ion_number = value
        self.update_last_interaction_filter()

    def update_last_interaction_filter(self):
        if self.packet_filter_mode == 'packet_nu':
            packet_filter = (self.last_line_interaction_angstrom > self.wavelength_start) & \
                      (self.last_line_interaction_angstrom < self.wavelength_end)
        elif self.packet_filter_mode == 'line_in_nu':
            line_in_nu = self.lines.wavelength.ix[self.last_line_interaction_in_id].values
            packet_filter = (line_in_nu > self.wavelength_start.to(u.angstrom).value) & \
                (line_in_nu < self.wavelength_end.to(u.angstrom).value)

                
        self.last_line_in = self.lines.ix[self.last_line_interaction_in_id[packet_filter]]
        self.last_line_out = self.lines.ix[self.last_line_interaction_out_id[packet_filter]]

        if self.atomic_number is not None:
            self.last_line_in = self.last_line_in[self.last_line_in.atomic_number == self.atomic_number]
            self.last_line_out = self.last_line_out[self.last_line_out.atomic_number == self.atomic_number]

        if self.ion_number is not None:
            self.last_line_in = self.last_line_in[self.last_line_in.ion_number == self.ion_number]
            self.last_line_out = self.last_line_out[self.last_line_out.ion_number == self.ion_number]


        last_line_in_count = self.last_line_in.wavelength.groupby(level=0).count()
        last_line_out_count = self.last_line_out.wavelength.groupby(level=0).count()

        self.last_line_in_table = self.lines[['wavelength', 'atomic_number', 'ion_number', 'level_number_lower',
                                              'level_number_upper']].ix[last_line_in_count.index]
        self.last_line_in_table['count'] = last_line_in_count
        self.last_line_in_table.sort('count', ascending=False, inplace=True)
        self.last_line_out_table = self.lines[['wavelength', 'atomic_number', 'ion_number', 'level_number_lower',
                                              'level_number_upper']].ix[last_line_out_count.index]
        self.last_line_out_table['count'] = last_line_out_count
        self.last_line_out_table.sort('count', ascending=False, inplace=True)


    def update_last_interaction_line_in_nu_filter(self):
        pass

    def plot_wave_in_out(self, fig, do_clf=True, plot_resonance=True):
        if do_clf:
            fig.clf()
        ax = fig.add_subplot(111)
        wave_in = self.last_line_list_in['wavelength']
        wave_out = self.last_line_list_out['wavelength']

        if plot_resonance:
            min_wave = np.min([wave_in.min(), wave_out.min()])
            max_wave = np.max([wave_in.max(), wave_out.max()])
            ax.plot([min_wave, max_wave], [min_wave, max_wave], 'b-')

        ax.plot(wave_in, wave_out, 'b.', picker=True)
        ax.set_xlabel('Last interaction Wave in')
        ax.set_ylabel('Last interaction Wave out')

        def onpick(event):
            print "-" * 80
            print "Line_in (%d/%d):\n%s" % (
                len(event.ind), self.current_no_packets, self.last_line_list_in.ix[event.ind])
            print "\n\n"
            print "Line_out (%d/%d):\n%s" % (
                len(event.ind), self.current_no_packets, self.last_line_list_in.ix[event.ind])
            print "^" * 80

        def onpress(event):
            pass

        fig.canvas.mpl_connect('pick_event', onpick)
        fig.canvas.mpl_connect('on_press', onpress)

class TARDISHistory(object):
    """
    Records the history of the model
    """


    def __init__(self, hdf5_fname, iterations=None):
        self.hdf5_fname = hdf5_fname
        if iterations is None:
            iterations = []
            hdf_store = pd.HDFStore(self.hdf5_fname, 'r')
            for key in hdf_store.keys():
                if key.split('/')[1] == 'atom_data':
                    continue
                iterations.append(int(re.match('model(\d+)', key.split('/')[1]).groups()[0]))

            self.iterations = np.sort(np.unique(iterations))
            hdf_store.close()
        else:
            self.iterations=iterations

        self.levels = None
        self.lines = None



    def load_atom_data(self):
        if self.levels is None or self.lines is None:
            hdf_store = pd.HDFStore(self.hdf5_fname, 'r')
            self.levels = hdf_store['atom_data/levels']
            self.lines = hdf_store['atom_data/lines']
            hdf_store.close()


    def load_t_inner(self, iterations=None):
        t_inners = []
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')

        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]

        for iter in iterations:
            t_inners.append(hdf_store['model%03d/configuration' %iter].ix['t_inner'])
        hdf_store.close()

        t_inners = np.array(t_inners)
        return t_inners

    def load_t_rads(self, iterations=None):
        t_rads_dict = {}
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')

        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]


        for iter in iterations:
            current_iter = 'iter%03d' % iter
            t_rads_dict[current_iter] = hdf_store['model%03d/t_rads' % iter]

        t_rads = pd.DataFrame(t_rads_dict)
        hdf_store.close()
        return t_rads

    def load_ws(self, iterations=None):
        ws_dict = {}
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')

        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]


        for iter in iterations:
            current_iter = 'iter%03d' % iter
            ws_dict[current_iter] = hdf_store['model%03d/ws' % iter]

        hdf_store.close()

        return pd.DataFrame(ws_dict)

    def load_level_populations(self, iterations=None):
        level_populations_dict = {}
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')
        is_scalar = False
        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            is_scalar = True
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]


        for iter in iterations:
            current_iter = 'iter%03d' % iter
            level_populations_dict[current_iter] = hdf_store['model%03d/level_populations' % iter]

        hdf_store.close()
        if is_scalar:
            return pd.DataFrame(level_populations_dict.values()[0])
        else:
            return pd.Panel(level_populations_dict)

    def load_jblues(self, iterations=None):
        jblues_dict = {}
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')
        is_scalar = False
        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            is_scalar = True
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]


        for iter in iterations:
            current_iter = 'iter%03d' % iter
            jblues_dict[current_iter] = hdf_store['model%03d/j_blues' % iter]

        hdf_store.close()
        if is_scalar:
            return pd.DataFrame(jblues_dict.values()[0])
        else:
            return pd.Panel(jblues_dict)


    def load_ion_populations(self, iterations=None):
        ion_populations_dict = {}
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')

        is_scalar = False
        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            is_scalar = True
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]


        for iter in iterations:
            current_iter = 'iter%03d' % iter
            ion_populations_dict[current_iter] = hdf_store['model%03d/ion_populations' % iter]

        hdf_store.close()
        if is_scalar:
            return pd.DataFrame(ion_populations_dict.values()[0])
        else:
            return pd.Panel(ion_populations_dict)

    def load_spectrum(self, iteration, spectrum_keyword='luminosity_density'):
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')

        spectrum = hdf_store['model%03d/%s' % (self.iterations[iteration], spectrum_keyword)]
        hdf_store.close()
        return spectrum

    def calculate_relative_lte_level_populations(self, species, iteration=-1):
        self.load_atom_data()
        t_rads = self.load_t_rads(iteration)
        beta_rads = 1 / (constants.k_B.cgs.value * t_rads.values[:,0])

        species_levels = self.levels.ix[species]

        relative_lte_level_populations = (species_levels.g.values[np.newaxis].T / float(species_levels.g.ix[0])) * \
                                         np.exp(-beta_rads * species_levels.energy.values[np.newaxis].T)

        return pd.DataFrame(relative_lte_level_populations, index=species_levels.index)

    def calculate_departure_coefficients(self, species, iteration=-1):
        self.load_atom_data()
        t_rads = self.load_t_rads(iteration)
        beta_rads = 1 / (constants.k_B.cgs.value * t_rads.values[:,0])

        species_levels = self.levels.ix[species]
        species_level_populations = self.load_level_populations(iteration).ix[species]
        departure_coefficient = ((species_level_populations.values * species_levels.g.ix[0]) /
                                 (species_level_populations.ix[0].values * species_levels.g.values[np.newaxis].T)) \
                                * np.exp(beta_rads * species_levels.energy.values[np.newaxis].T)

        return pd.DataFrame(departure_coefficient, index=species_levels.index)

    def get_last_line_interaction(self, iteration=-1):
        iteration = self.iterations[iteration]
        self.load_atom_data()

        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')
        model_string = 'model'+('%03d' % iteration) +  '/%s'
        last_line_interaction_in_id = hdf_store[model_string % 'last_line_interaction_in_id'].values
        last_line_interaction_out_id = hdf_store[model_string % 'last_line_interaction_out_id'].values
        last_line_interaction_shell_id = hdf_store[model_string % 'last_line_interaction_shell_id'].values
        try:
            montecarlo_nu = hdf_store[model_string % 'montecarlo_nus_path'].values
        except KeyError:
            montecarlo_nu = hdf_store[model_string % 'montecarlo_nus'].values
        hdf_store.close()
        return LastLineInteraction(last_line_interaction_in_id, last_line_interaction_out_id, last_line_interaction_shell_id,
                            montecarlo_nu, self.lines)
