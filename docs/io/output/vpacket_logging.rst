**********************
Virtual Packet Logging
**********************

If ``virtual_packet_logging`` is set to ``True`` in either the :ref:`spectrum configuration <spectrum-config>` or as 
an :ref:`argument to the run_tardis function <optional-input>`, then TARDIS will log information about the
simulation's :ref:`virtual packets <virtual_packets>`.

After running the simulation, the following information can be retrieved:

.. note::
    The following tables list the attributes of the simulation in which the information would be stored. For
    example, if you call your TARDIS simulation ``sim``, you would access the virtual packet frequencies by running
    ``sim.runner.virt_packet_nus``.


.. list-table::
    :header-rows: 1
 
    * - Attribute of Simulation
      - Description
    * - ``virt_packet_logging``
      - Boolean that shows if virtual packet logging is turned on
    * - ``runner.virt_packet_nus``
      - List of virtual packet frequencies
    * - ``runner.virt_packet_energies``
      - List of virtual packet energies
    * - ``runner.virt_packet_last_interaction_type``
      - Type of interaction that caused the virtual packet to be spawned
    * - ``runner.virt_packet_last_interaction_in_nu``
      - The frequency the virtual packet was spawned at
    * - ``runner.virt_packet_last_line_interaction_in_id``
      - If the last interaction was a line interaction, the line_interaction_in_id for that interaction (see :doc:`physical_quantities`)
    * - ``runner.virt_packet_last_line_interaction_out_id``
      - If the last interaction was a line interaction, the line_interaction_out_id for that interaction (see :doc:`physical_quantities`)