#codes to for analyse the model.

from astropy import units as u


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