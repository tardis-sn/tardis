*************
Configuration
*************

The TARDIS configuration consists of multiple sections that pertain to certain parts of the code. We will use the
schemas to show what options are available. Our schema mark-up defines names in bold-fat as required.
can be seen here:

Base Schema
===========

.. jsonschema:: schemas/base.yml


The base schema links to several other configuration parts

Sub Schemas
===========

.. jsonschema:: schemas/supernova.yml

.. jsonschema:: schemas/model.yml

.. jsonschema:: schemas/plasma.yml

.. jsonschema:: schemas/montecarlo.yml

.. jsonschema:: schemas/spectrum.yml