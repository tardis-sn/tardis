*******************************
Hierarchical Data Format (HDF5)
*******************************


Overview
========


What is HDF5?
-------------

`HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ is a high-performance data software library and file format to manage, process, and store your heterogeneous data. HDF5 is built for fast I/O processing and storage.

In TARDIS, it is used to store the data from simulations and other actions.


HDF5 structure
--------------

All this data is stored in the following structure:

- **File:** a contiguous string of bytes in a computer store (memory, disk, etc.), and the bytes represent zero or more objects of the model.
- **Group:** a collection of objects (including groups).
- **Dataset:** a multidimensional array of data elements with attributes and other metadata.
- **Dataspace:** a description of the dimensions of a multidimensional array.
- **Datatype:** a description of a specific class of data element, including its storage layout as a pattern of bits.
- **Attribute:** a named data value associated with a group, dataset, or named datatype.
- **Property List:** a collection of parameters (some permanent and some transient) controlling options in the library.
- **Link:** the way objects are connected

This structure is stored in a file with the extension “.h5”.


Libraries
---------

TARDIS uses `these libraries <https://www.hdfgroup.org/downloads/hdf5>`_ in Python, but they could be installed on the computer to work with the H5 files. The advantage of installing the libraries on the computer is that it allows you to handle the HDF5 with command-line instructions quickly and powerfully.


HDFView
-------

`HDFView is a visual tool <https://www.hdfgroup.org/downloads/hdfview/>`_ written in Java for browsing and editing HDF (HDF5 and HDF4) files. Using HDFView, you can:

- View a file hierarchy in a tree structure.
- Create new files; add or delete groups and datasets.
- View and modify the content of a dataset.
- Add, delete, and modify attributes.


Extract data
============


Extract data using the HDF5 libraries
-------------------------------------

In TARDIS the generated data is valuable because it has been reviewed and validated. For this reason, TARDIS tests against the generated data.

Then the question is how to extract data from the actual “.h5” files. It could be easily pulled and followed by the next commands:

Copy the specific group of data in a new “.h5” file.

.. code-block:: shell

    h5copy \
     --input unit_test_data.h5 \
     --source "/test_transport_simple/transport_state" \
     --output TestTransportSimple.h5\
     --destination "/transport_state"

Compare and validate that the original and copy data are the same.

.. code-block:: shell

    h5diff \
     -r \
     unit_test_data.h5 \
     TestTransportSimple.h5 \
     /test_transport_simple/transport_state/j_blue_estimator \
     /transport_state/j_blue_estimator

It produces the following output, which confirms that the data is the same:

.. code-block:: text

    attribute: <CLASS of </test_transport_simple/transport_state/j_blue_estimator>> and <CLASS of </transport_state/j_blue_estimator>>
    0 differences found
    attribute: <TITLE of </test_transport_simple/transport_state/j_blue_estimator>> and <TITLE of </transport_state/j_blue_estimator>>
    0 differences found
    …
    attribute: <VERSION of </test_transport_simple/transport_state/j_blue_estimator/block0_values>> and <VERSION of </transport_state/j_blue_estimator/block0_values>>
    0 differences found
    attribute: <transposed of </test_transport_simple/transport_state/j_blue_estimator/block0_values>> and <transposed of </transport_state/j_blue_estimator/block0_values>>
    0 differences found
