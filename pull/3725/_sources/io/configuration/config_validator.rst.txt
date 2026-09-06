***********************
Configuration Validator
***********************

The default config validator takes a user configuration and a default configuration and creates a consistent and valid configuration for TARDIS based on the constraints given in the default configuration. Both input data are normally given as a YAML dictionary with a consistent hierarchical structure, i.e. for every item in the user configuration there has to be a declaration in the default configuration at the same hierarchical level. This declaration can be either an unspecific empty level declaration like:

- Main_level:
	- Second_level:
		- Third_level:
			…

or a declaration of a configuration item like:

- item:
    - property_type: int
    - default: 1
    - mandatory: True
    - help:  ‘This is a doc string.'

This contains always  the keywords ``help``, ``default``, ``mandatory``, and ``property_type``. The keyword help is  a doc-string which describes the corresponding item. "Default" specifies the default value which is used in case that no value for this item is specified in the corresponding user configuration item. If the keyword mandatory is ``True``, the item has to be specified in the user configuration. The keyword ``property_type`` is used to specify the type of the item. At the moment, the config parser knows the following types:

.. ::

- **Int:** The property type int is for integer like config items.
- **Float:** The property type float is for float like config items.
- **String:** The property type string is for string like config items.
- **Quantity:** The property type quantity is for physical quantities with units given as string. The string contains value and unit separated by a whitespace E.g. 2 cm.
- **Range:** The property type range specifies a range via start and end.

.. note:: ``abs(start - end ) > 0``.

- **Quantity_range:** Like property type range but with quantities as start and stop. The consistency of the units is checked.

Additionally to the four standard keywords the types integer, float, and quantity can have the keywords ``allowed_value`` and ``allowed_type``. ``allowed_value`` specifies the allowed values in a list, whereas ``allowed_type`` specifies a range of allowed values, such as “x>10”.

Container
^^^^^^^^^

For more complex configurations with dependencies, you can use the containers that allow branching in the configuration. A container is declared in the default configuration file by setting the  ``property_type`` to container property and specifying the properties of the container with keyword type. The ``property_type`` of this section is container-declaration that allows you to specify the possible container items with the keyword container. For every specified container item, the code expects the declaration of all sub items. The keywords for this are “_“ + “name of the container item”.
If the type declaration for this container is finished, you can specify all container items like normal items. Here is an example for a container configuration with two branches:

.. source: yaml

- container_example:
        - property_type: container-property
        - type:
            - property_type: container-declaration
            - containers: ['one', 'two', 'three']
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
                - default: [1 m,10 cm] *#[Start,End]*
                - mandatory: False
                - help:  Like property type range but with quantities as start and stop. The consistency of the units is checked.

How to use
^^^^^^^^^^
                
To use the default parser create a new config object form the class ConfigurationValidator either from a dictionaries or from YAML files:

::

- My_config = ConfigurationValidator(default configuration dictionary, user configuration dictionary)

or

::

- My_config = ConfigurationValidator.from_yaml (default configuration file, user configuration file).

To access the configuration for TARDIS, use the method ``get_config``.

