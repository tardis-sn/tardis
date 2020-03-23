#### ONLY FOR QUICKLY TESTING IMPROVEMENTS WHILE WORKING
#### TESTED ON CHANGES IN PROGRESS. IF RUN, MAY GIVE ERROR
#### TO BE EXCLUDED FROM FINAL MERGE COMMIT

from tardis.montecarlo.montecarlo_numba import r_packet
import time

start = time.time()
print(r_packet.calculate_distance_boundary(1, 2, 3, 4))
t = time.time()-start
print("Time for execution: {}".format(t))

# rpack = r_packet.RPacket(1, 2, 3, 4, 5)
# start = time.time()
# print(r_packet.calculate_distance_line(rpack, 5, 2, 3))
# t = time.time()-start
# print("Time for execution: {}".format(t))