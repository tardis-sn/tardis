import tardis.util.base as base
import time

start = time.time()
val2 = base.intensity_black_body(1,5)
time_diff = time.time() - start

print("Exec Time: {}".format(time_diff))