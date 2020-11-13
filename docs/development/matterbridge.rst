************
Matterbridge
************

========
Overview
========

We use Matterbridge to connect channels across different messaging platforms
such as Slack, Mattermost and Gitter. Matterbridge provides binaries for many
operating systems. These files are intended to be used in combination with
TOML configuration files::
  $ ./matterbridge -conf config-tardis-matterbridge.toml

The TOML file includes all the configuration needed to connect as many
messaging platforms as you want, including tokens and passwords.

Currently, we keep a service running on `OpenSupernova.org`_ server to keep
running Matterbridge as a daemon. Our configuration files are stored in the
private repository ``tardis-matterbridge``. These files include our custom
TOML configuration file and the files needed to set up the Linux service.

Since this server runs Ubuntu 14.04 we use an `Upstart script`_ instead a
`Systemd service`_ file. But a *Systemd* file is also included in
``tardis-matterbridge`` repository, just in case.


=====
Setup
=====

Use the TOML file in ``tardis-matterbridge`` as an example. You can set up
as many gateways as you want.

-----
Slack
-----

Follow the `steps for Slack in the Matterbridge wiki`_ and read carefully
`this comment`_. Do not add any scopes manually. Copy the token to use later
in the TOML file.


------
Gitter
------

Create a new GitHub dedicated account to use as a bot, refer to `this link`_
and copy the token.


----------------------------------------
Setting up the server for the first time
----------------------------------------

Follow these steps to set up the server for the first time:

1. *ssh* to `OpenSupernova.org`_
2. Download the Matterbridge binary for Linux from the `releases section`_
3. Make the file executable and rename it to ``matterbridge``
4. Copy ``matterbridge`` executable to ``/usr/local/bin``
5. Clone ``tardis-matterbridge`` repository
6. Copy ``config-tardis-matterbridge.toml`` to ``/usr/local/etc/matterbridge/``
7. Copy ``matterbridge.conf`` to ``/etc/init/``
8. Run ``sudo service matterbridge start``
9. Test your gateways


-------------------
Adding new gateways
-------------------

After updating the TOML file, follow these steps:

1. *ssh* to `OpenSupernova.org`_
2. Copy your new ``config-tardis-matterbridge.toml`` to ``/usr/local/etc/matterbridge/``
3. Run ``sudo service matterbridge restart``
4. Test your gateways
5. If everything is ok, make a pull request to ``tardis-matterbridge`` with your new TOML file

.. include:: matterbridge.inc
