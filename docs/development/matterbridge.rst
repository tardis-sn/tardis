*******************
Matterbridge Server
*******************

========
Overview
========

We use Matterbridge to connect channels across different messaging platforms
such as *Slack*, *Mattermost* and *Gitter*. Matterbridge provides binaries 
for  many operating systems. The ``matterbridge`` binary is intended to be used

in combination with a :term:`TOML` configuration file::

  $ ./matterbridge -conf config-tardis-matterbridge.toml

The TOML file includes all the parameters required to connect as many
messaging platforms as you want, like *tokens* and *passwords*. Once the
application is running, messages can be shared between the connected rooms.

Currently, we keep a service running on the `OpenSupernova.org`_ server to 
run Matterbridge as a daemon. Our configuration files are stored in a `private 
GitHub repository`_, including our custom TOML and the files needed to set up the
Linux service.

Since this server runs Ubuntu 14.04 we use an `Upstart script`_ instead of a
`Systemd service`_ file, but a Systemd file is also included in our repository.
A *cron* job restarts the service periodically to prevent disconnections.


=============
Configuration
=============

Use the TOML file in ``tardis-matterbridge`` as an example. You can set up
as many gateways as you want!


-----
Slack
-----

Follow the `steps for Slack in the Matterbridge wiki`_ and read carefully
`this comment`_ (do not add any scopes manually!). Copy the token to use
later in the TOML file.


------
Gitter
------

Create a new GitHub dedicated account to use as a bot, refer to `this link`_
and copy the token to use later in the TOML file.


----------
Mattermost
----------

Follow the `steps for Mattermost in the Matterbridge wiki`_.


================
First-time setup
================

Follow these steps to set up the server:

1. *ssh* to `OpenSupernova.org`_ server.
2. Download the Matterbridge binary for Linux from the `releases section`_.
3. Make the file executable and rename it to ``matterbridge``.
4. Copy ``matterbridge`` executable to ``/usr/local/bin``.
5. Clone ``tardis-matterbridge`` repository in your ``$HOME``.
6. Copy ``config-tardis-matterbridge.toml`` to ``/usr/local/etc/matterbridge/``.
7. Copy ``matterbridge.conf`` to ``/etc/init/``.
8. Run ``sudo service matterbridge start``.
9. Test your gateways.


===========================
Update server configuration
===========================

After updating the TOML file, follow these steps:

1. *ssh* to `OpenSupernova.org`_ server.
2. Copy your new ``config-tardis-matterbridge.toml`` to ``/usr/local/etc/matterbridge/``.
3. Run ``sudo service matterbridge restart``.
4. Test your gateways.
5. If everything is ok, make a pull request to ``tardis-matterbridge`` with your new TOML file.


========
Cron job
========

To edit the cron job, *ssh* to `OpenSupernova.org`_ and run ``sudo crontab -e``.


.. include:: matterbridge.inc
