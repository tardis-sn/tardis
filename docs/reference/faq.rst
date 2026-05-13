.. _faq:

**************************
Frequently Asked Questions
**************************

**This reference section provides quick answers to common questions about TARDIS.**

This FAQ serves as a reference for finding solutions to frequently encountered issues. For step-by-step learning, see our :doc:`tutorials`. For specific problem-solving guides, check our :doc:`how_to_guides`.

Overview
--------

Here you'll find answers to common questions organized by topic.

- :ref:`faq-usage`

     - :ref:`faq-usage-memory`


.. _faq-usage:

Usage
-----

.. _faq-usage-memory:

My simulation seems to consume excessive amounts of memory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


High memory usage can have many reasons. Both your model
and certain simulation settings can increase the memory
usage significantly

1. Enabling ``track_rpacket: true`` will take up substantial
   amounts of memory, in particular for higher packet counts.
   Consider turning this feature off or only use it with a
   small amount of :term:`packets`.
2. Both the number of shells and the number of :term:`packets`
   increase the memory requirements. Consider using only
   as many shells/ :term:`packets` as are required to converge
   on a result.
