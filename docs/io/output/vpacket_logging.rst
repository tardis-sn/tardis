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
      - Type
      - Description
    * - ``virt_packet_logging``
      - Boolean
      - Shows if virtual packet logging is turned on
    * - ``runner.virt_packet_nus``
      - Numpy array
      - Virtual packet frequencies
    * - ``runner.virt_packet_energies``
      - Numpy array
      - Virtual packet energies
    * - ``runner.virt_packet_initial_mus``
      - Numpy array
      - Propagation directions that virtual packets are launched at
    * - ``runner.virt_packet_initial_rs``
      - Numpy array
      - Radii that virtual packets are launched at
    * - ``runner.virt_packet_last_interaction_type``
      - Numpy array
      - | Type of interaction that caused the virtual packets to be spawned
        | (enum, see :doc:`physical_quantities`)
    * - ``runner.virt_packet_last_interaction_in_nu``
      - Numpy array
      - Frequencies of the r-packets which spawned the virtual packet
    * - ``runner.virt_packet_last_line_interaction_in_id``
      - Numpy array
      - | If the last interaction was a line interaction, the
        | line_interaction_in_id for that interaction 
        | (see :doc:`physical_quantities`)
    * - ``runner.virt_packet_last_line_interaction_out_id``
      - Numpy array
      - | If the last interaction was a line interaction, the
        | line_interaction_out_id for that interaction 
        | (see :doc:`physical_quantities`)