Matterbridge
------------

.. _how-to-guide-edit-the-cron-job:

How-To Guide: Edit the Cron Job
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SSH to OpenSupernova.org and run:

.. code-block:: shell

   sudo crontab -e


.. _explanation-matterbridge:

Explanation: Matterbridge
~~~~~~~~~~~~~~~~~~~~~~~~~

Matterbridge connects messaging channels across platforms such as Slack,
Mattermost, and Gitter. The ``matterbridge`` binary is used with a TOML
configuration file:

.. code-block:: shell

   ./matterbridge -conf config-tardis-matterbridge.toml


The TOML file contains parameters required to connect rooms, including tokens
and passwords. When the application runs, messages can be shared between the
connected rooms.

TARDIS keeps a service running on the OpenSupernova.org server to run
Matterbridge as a daemon. Configuration files are stored in a private GitHub
repository, including the custom TOML and Linux service files.

The server runs Ubuntu 14.04, so TARDIS uses an Upstart script instead of a
Systemd service file. A Systemd file is also included in the repository. A cron
job restarts the service periodically to prevent disconnections.

.. _reference-matterbridge-reference:

Reference: Matterbridge Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Matterbridge command:

.. code-block:: shell

   ./matterbridge -conf config-tardis-matterbridge.toml


Matterbridge configuration:

- Use the TOML file in ``tardis-matterbridge`` as an example.
- Configure as many gateways as needed.
- For Slack, follow the Matterbridge wiki Slack setup steps and read the linked
  comment warning not to add scopes manually.
- For Gitter, create a dedicated GitHub bot account and copy the token.
- For Mattermost, follow the Matterbridge wiki Mattermost setup steps.

Important paths:

- Matterbridge executable: ``/usr/local/bin/matterbridge``
- Matterbridge TOML config:
  ``/usr/local/etc/matterbridge/config-tardis-matterbridge.toml``
- Upstart service config: ``/etc/init/matterbridge.conf``

Important repositories and services:

- OpenSupernova.org server: http://opensupernova.org
- Private Matterbridge repository:
  https://github.com/tardis-sn/tardis-matterbridge
- Matterbridge releases: https://github.com/42wim/matterbridge/releases/latest

.. _how-to-guide-set-up-matterbridge:

How-To Guide: Set Up Matterbridge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


1. SSH to the OpenSupernova.org server.
2. Download the Matterbridge binary for Linux from the releases section.
3. Make the file executable and rename it to ``matterbridge``.
4. Copy the ``matterbridge`` executable to ``/usr/local/bin``.
5. Clone the ``tardis-matterbridge`` repository in ``$HOME``.
6. Copy ``config-tardis-matterbridge.toml`` to
   ``/usr/local/etc/matterbridge/``.
7. Copy ``matterbridge.conf`` to ``/etc/init/``.
8. Start the service:

   .. code-block:: shell

      sudo service matterbridge start


9. Test the gateways.

.. _how-to-guide-update-server-configuration:

How-To Guide: Update Server Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After updating the TOML file:

1. SSH to the OpenSupernova.org server.
2. Copy the new ``config-tardis-matterbridge.toml`` to
   ``/usr/local/etc/matterbridge/``.
3. Restart the service:

   .. code-block:: shell

      sudo service matterbridge restart


4. Test the gateways.
5. If everything works, open a pull request to ``tardis-matterbridge`` with the
   new TOML file.

Matterbridge links:

- OpenSupernova.org: http://opensupernova.org
- Private Matterbridge repository:
  https://github.com/tardis-sn/tardis-matterbridge
- Upstart script:
  https://www.digitalocean.com/community/tutorials/the-upstart-event-system-what-it-is-and-how-to-use-it
- Systemd service:
  https://freedesktop.org/software/systemd/man/systemd.service.html
- Slack setup:
  https://github.com/42wim/matterbridge/wiki/Slack-bot-setup
- Mattermost setup:
  https://github.com/42wim/matterbridge/wiki/Section-Mattermost-%28basic%29
- Slack scopes comment:
  https://github.com/42wim/matterbridge/issues/964#issuecomment-612721850
- Gitter authentication:
  https://developer.gitter.im/docs/authentication
- Matterbridge releases:
  https://github.com/42wim/matterbridge/releases/latest


