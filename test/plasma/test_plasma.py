from tardis import atomic, plasma
import sqlite3
import numpy as np

conn = sqlite3.connect('/Users/wkerzend/tmp/tardis_test/test.db3', detect_types=sqlite3.PARSE_DECLTYPES)

atomic_model = atomic.CombinedAtomicModel.from_db(conn)
abundances = np.zeros(30)
abundances[19]=1.

sn_plasma = plasma.NebularPlasma.from_model({'ca':1.}, 1e-17, atomic_model)
conn.close()