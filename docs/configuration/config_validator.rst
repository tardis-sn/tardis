***********************
Configuration Validator
***********************

.. warning::
  The configuration validator along with the syntax described here is subject to change soon. This document serves as a reference for the developers of TARDIS.

The default configuration validator takes a configuration definition and a user configuration and creates a consistent and valid configuration for tardis based on the constraints given in the configuration definition.

Both input data are normally given as a yaml dictionary with a consistent hierarchical structure, i.e. for every item in the user configuration there has to be a declaration in the default configuration  at the same hierarchical level. This declaration can be either an unspecific empty level declaration like:

.. code-block:: yaml

  - first_level:
    - second_level:
      - third_level:
        …

or a declaration of a configuration property like:

.. code-block:: yaml

  - my_property:
    - property_type: int
    - default: 1
    - mandatory: True
    - help: This is a doc string.
        
Property definitions must always contain the following keywords:

property_type
  The expected type of the property. See `Property Types`_
default
  The default value which the property gets when no value is provided in the user configuration.
mandatory
  If set to True, the property is required to be specified in the user configuration.
help
  A doc-string which describes the property.


Property Types
^^^^^^^^^^^^^^^

At the moment, the default configuration makes use of the following property types:

* *int*: For integer properties.
* *float*: For floating point properties.
* *bool*: For boolean (True/False) properties.
* *string*: For text properties.
* *quantity*: For physical quantities expressed as a value and a unit separated by a whitespace (e.g.: 2 cm)
* *quantity_range_sampled*: For specifying sampled ranges of quantities from ``start`` to ``end``, with ``num`` samples. ``start`` and ``end`` should be quantities, while ``num`` should be an integer. The consistency of the units is checked.
* *list*: A list of values.
* *container-property* / *container-declaration*: For creating custom complex virtual types. See `The container property`_

.. note::
  The types ``int``, ``float``, and ``quantity`` except from the standard keywords can also have the keywords ``allowed_value`` and ``allowed_type``. ``allowed_value`` specifies the allowed values of the property in a list, whereas ``allowed_type`` specifies a range of allowed values for the property (e.g. ``x>10``).

The configuration validator also supports a few more types, which are not currently used by TARDIS. These are:

* *range*: For specifying ranges from ``start`` to ``end``. (Note: :code:`abs(start - end) > 0`)
* *range_sampled*: For specifying sampled ranges from ``start`` to ``end``, with ``num`` samples. (Note: :code:`abs(start - end) > 0`)
* *quantity_range*: For specifying ranges of quantities from ``start`` to ``end``. Both ``start`` and ``end`` should be quantities and the consistency of the units is checked.
* *abundance_set*
* *legacy-abundances*



The container property
^^^^^^^^^^^^^^^^^^^^^^

For more complex configurations with dependencies, the container type can be used which allows creating a complex property consisting of other basic ones.

A container is declared in the configuration definition file by setting the ``property_type`` of the property which should be a container to ``container-property``. Then the properties of the container can be specified inside the keyword ``type``.

.. note::
The ``property_type`` of the ``type`` section of a container property should be ``container-declaration``.
  
The ``containers`` keyword in a ``container-declaration`` specifies a list of all the allowed values a user can give as the type of a property in the user configuration. For every virtual type in the ``containers`` list there must be a declaration of what sub-properties this type allows or requires.

The syntax for declaring the sub-properties is a keyword for each defined type formed as: ``_type_name:`` followed by a list of the **mandatory** sub-properties that should be provided with this type, and, optionally, another keyword formed as ``+type_name:`` followed by a list of the **optional** sub-properties that can be provided with this type.

After specifying the sub-properties that should or could be provided with each type, the sub-properties can be defined like usual properties.

Following is an example for a container property with two possible types, ``one`` and ``two``:

.. code-block:: yaml

  - container_example:
    - property_type: container-property
    - type:
      - property_type: container-declaration
      - containers: ['one', 'two']
      - _one: ['one_one', 'one_two']
      - _two: ['two_one']

    - one_one:
      - property_type: string
      - default: 'This is a container item'
      - mandatory: False
      - help: This is a container item from the container one.
    
    - one_two:
      - sub_one_two_one:
        - property_type: string
        - default: 'This is a container item'
        - mandatory: False
        - help: This is a container item from the container one.
      - sub_one_two_two:
        - property_type: string
        - default: 'This is a container item'
        - mandatory: False
        - help: This is a container item from the container one.
    
    - two_one:
      - quantity_range:
        - property_type: quantity_range
        - default: [1 m,10 cm] #[Start,End]
        - mandatory: False
        - help:  Like property type range but with quantities as start and stop. The consistency of the units is checked.
