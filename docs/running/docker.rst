Running TARDIS in a docker container
====================================
Before getting started, make sure you have docker installed from the
`official website <https://docs.docker.com/install/>`_

**To setup and run the container (under development) follow these steps:**

**1**. Pull the image from the relevant repository
The TARDIS docker container currently resides in `dockerhub <https://hub.docker.com/>`_.
Pull it from the repository
.. code-block:: shell

    docker pull not4win/tardis:1.0

Using with command-line
-----------------------
After pulling the image, run it on interactive mode with following commands.
.. code-block:: shell

    docker run -it not4win/tardis:1.0

Depending on your requirements, ``build`` or set the project in ``develop`` mode.
.. code-block:: shell

    python setup.py develop

You are good to go!.

Using with jupyter notebook
---------------------------
Steps are little bit different from the above one as there is a requirement of port forwarding from
docker to localhost. ``-p`` argument is used for port forwarding.
.. code-block:: shell

    docker run -it -p 8999:8999 not4win/tardis:1.0

Here ``8999`` is a port number needed for port forwarding while using jupyter notebook.

Depending on your requirements, ``build`` or set the project in ``develop`` mode.
.. code-block:: shell

    python setup.py develop

Run the jupyter notebook.
.. code-block:: shell

    jupyter notebook --port=8999 --ip=0.0.0.0 --allow-root

You are good to go! For an exercise, just run the `quickstart-code <../quickstart/quickstart.ipynb>`_



