**********************
Virtual Packet Logging
**********************

If ``virtual_packet_logging`` is set to ``True`` in either the :ref:`spectrum configuration <spectrum-config>` or as 
an :ref:`argument to the run_tardis function <optional-input>`, then TARDIS will log information about the
simulation's :ref:`virtual packets <virtual_packets>`.

The parent-interaction fields are populated from the real-packet histories
used by virtual-packet spectrum post-processing. Virtual packets spawned at the
inner boundary have ``last_interaction_type`` ``NO_INTERACTION``, ``before_nu``
0, ``radius`` equal to the inner boundary radius, and -1 for ``shell_id`` and
the line ids.

After running the simulation, the logged virtual packets are available from the
spectrum solver as a DataFrame with one row per escaping virtual packet. For a
TARDIS simulation ``sim``, run ``sim.spectrum_solver.vpacket_tracker_to_df()``;
workflows provide the same method on ``workflow.spectrum_solver``.
The interaction columns use the same names and types as the real-packet
``transport_state.tracker_last_interaction_df``. Each row
describes the real-packet event that spawned the virtual packet.

.. list-table::
    :header-rows: 1

    * - Column
      - Type
      - Description
    * - ``last_interaction_type``
      - Categorical
      - | Type of the real-packet event that spawned the virtual packet;
        | ``NO_INTERACTION`` for virtual packets spawned at the inner boundary
    * - ``radius``
      - Float
      - Radius of the spawning real-packet event [cm]
    * - ``shell_id``
      - Integer
      - Shell of the spawning real-packet event
    * - ``before_nu``
      - Float
      - Frequency of the real packet before the spawning event [Hz]
    * - ``line_absorb_id``
      - Integer
      - | Absorbing line of the spawning event, or -1 if it was not a line
        | interaction (see :doc:`../../how_to/output/how_to_physical_quantities`)
    * - ``line_emit_id``
      - Integer
      - | Emitting line of the spawning event, or -1 if it was not a line
        | interaction (see :doc:`../../how_to/output/how_to_physical_quantities`)
    * - ``nus``
      - Float
      - Escaping virtual packet frequency [Hz]
    * - ``energies``
      - Float
      - Escaping virtual packet energy [erg]
    * - ``initial_mus``
      - Float
      - Propagation direction the virtual packet was launched at
    * - ``initial_rs``
      - Float
      - Radius the virtual packet was launched at [cm]

When a simulation is saved to HDF, this DataFrame is stored at
``/simulation/spectrum_solver/vpacket_tracker_df``, with
``last_interaction_type`` stored as ``InteractionType`` integer codes.
