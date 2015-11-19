import pprint
import numpy as np


class BoundFreeException(Exception):
    pass


class IncompletePhotoionizationDataError(BoundFreeException):
    def __init__(self, message=None, needed_data=None, provided_data=None, list_data_mismatch=False):
        if message == None:
            message = 'Not all needed photoionization data was provided.'
        super(IncompletePhotoionizationDataError, self).__init__(message)
        self.needed_photoion_data = needed_data
        self.provided_photoion_data = provided_data
        if list_data_mismatch == True:
            self._list_needed_and_provided_data()

    def _list_needed_and_provided_data(self):
        opt = np.get_printoptions()
        np.set_printoptions(threshold='nan')
        print 'Provided photoionization data for levels (i,j,k):'
        pprint.pprint(self.provided_photoion_data)
        print 'Needed photoionization data for levels (i,j,k):'
        pprint.pprint(self.needed_photoion_data)
        np.set_printoptions(**opt)