#codes to for analyse the model.

import re
import os

from astropy import units as u, constants
import numpy as np
import pandas as pd


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

    # print "Lines in:\n------\n %s" % line_in_table
    # print "Lines out:\n------\n %s" % line_out_table

    return last_line_in_ids, last_line_out_ids


class LastLineInteraction(object):

    def __init__(self, model, packet_filter_mode='packet_nu'):
        self.model = model
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
            packet_filter = (self.model.last_line_interaction_angstrom > self.wavelength_start) & \
                      (self.model.last_line_interaction_angstrom < self.wavelength_end)
        elif self.packet_filter_mode == 'line_in_nu':
            line_in_nu = self.model.atom_data.lines.wavelength.ix[self.model.last_line_interaction_in_id].values
            packet_filter = (line_in_nu > self.wavelength_start.to(u.angstrom).value) & \
                (line_in_nu < self.wavelength_end.to(u.angstrom).value)

        self.last_line_in = self.model.atom_data.lines.ix[self.model.last_line_interaction_in_id[packet_filter]]
        self.last_line_out = self.model.atom_data.lines.ix[self.model.last_line_interaction_out_id[packet_filter]]



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


    def __init__(self, hdf5_fname, history_dir='.', iterations=None):
        self.hdf5_fname = os.path.join(history_dir, hdf5_fname)
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

        if iterations is None:
            iterations = self.iterations
        elif np.isscalar(iterations):
            iterations = [self.iterations[iterations]]
        else:
            iterations = self.iterations[iterations]


        for iter in iterations:
            current_iter = 'iter%03d' % iter
            level_populations_dict[current_iter] = hdf_store['model%03d/level_populations' % iter]

        hdf_store.close()

        return pd.Panel(level_populations_dict)

    def calculate_departure_coefficients(self, iteration=-1):
        iteration = self.iterations[iteration]

        t_rads = self.load_t_rads(iteration)
        beta_rads = 1 / (constants.k_B.cgs.value * t_rads.values)


    def get_spectrum(self, iteration, spectrum_keyword='luminosity_density'):
        hdf_store = pd.HDFStore(self.hdf5_fname, 'r')

        spectrum = hdf_store['model%03d/%s' % (self.iterations[iteration], spectrum_keyword)]
        hdf_store.close()
        return spectrum

    def plot_convergence(self, fig):

        ax = fig.add_subplot(111)

        ax.plot(self.t_inner)
        ax.set_xlabel('iterations')
        ax.set_ylabel('t_inner [K]')


def get_departure_coefficient(fname, iteration, species, model_dir='.'):
    fname = os.path.join(model_dir, fname)
    level_populations = pd.HDFStore(fname,'r')['model%03d/level_populations' % iteration].ix[species]
    beta_rad = 1/(constants.k_B.cgs.value * pd.HDFStore(fname, 'r')['model%03d/t_rads' % iteration].values)
    levels = pd.HDFStore(fname, 'r')['atom_data/levels'].ix[species]

    departure_coefficient = ((level_populations.values * levels.g.ix[0]) / (level_populations.ix[0].values * levels.g.values[np.newaxis].T)) * np.exp(beta_rad * levels.energy.values[np.newaxis].T)
    t_rad = pd.HDFStore(fname, 'r')['model%03d/t_rads' % iteration].values
    w = pd.HDFStore(fname, 'r')['model%03d/ws' % iteration].values
    return t_rad, w, levels, departure_coefficient