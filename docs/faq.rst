.. _installation:

**************************
Frequently Asked Questions
**************************

Usage
-----

.. _faq-usage:

My simulation seems to consume excessive amounts of memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

High memory usage can have many reasons. Both your model
and certain simulation settings can increase the memory
usage significantly

1. Enabling ``track_rpacket: true`` will take up substantial
   amounts of memory, in particular for higher packet counts.
   Consider turning this feature off or only use it with a
   small amount of packets.
2. Both the number of shells and the number of packets
   increase the memory requirements. Consider using only
   as many shells/ packets as are required to converge
   on a result.
