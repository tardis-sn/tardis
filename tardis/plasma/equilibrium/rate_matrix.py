import pandas as pd


class RateMatrix:
    def __init__(self, rates: list):
        self.rates = rates

    def solve(self):
        rate_matrix = pd.DataFrame()
        for rate in self.rates:
            rate_df = rate.solve()
            rate_matrix += rate_df

        # get unique species (atomic, ion numbers) and group by them
        # fill NAN values with 0
        """ rate_matrices
        for species_id, rates in rates.groupby(level=(‘atomic_number’, ‘ion_number’)):
            rate_matrices.loc[species_id] = make_rate_matrix_from_rates """
        # make first n = 1 to normalize (see plasma docs and L237 partition_functions.py)

        return rate_matrix
