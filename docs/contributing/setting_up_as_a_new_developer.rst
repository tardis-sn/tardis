.. _new-setting-up-as-a-new-developer:

Setting Up As A New Developer
=============================

Use this page when you are preparing a local checkout for TARDIS development.
It gathers the one-time setup steps that used to be spread across the Git,
code-style, and codebase-structure pages.

.. _new-how-to-guide-set-up-github-access:

Set Up GitHub Access
--------------------

1. Create a GitHub account at https://github.com if needed.
2. Configure write access with SSH keys.
3. Use GitHub Help's SSH key instructions:
   https://help.github.com/articles/generating-ssh-keys

.. _new-how-to-guide-fork-and-clone-tardis:

Fork And Clone TARDIS
---------------------

You only need to fork once for each package you want to contribute to.

1. Log into GitHub.
2. Go to the TARDIS GitHub home page.
3. Click the fork button.
4. Clone your fork:

   .. code-block:: shell

      git clone git@github.com:your-user-name/tardis.git
      cd tardis

5. Check branches:

   .. code-block:: shell

      git branch -a

6. Check remotes:

   .. code-block:: shell

      git remote -v

At this point, ``origin`` should point to your GitHub fork.

.. _new-how-to-guide-connect-your-fork-to-upstream:

Connect Your Fork To Upstream
-----------------------------

Add the main TARDIS repository as ``upstream``:

.. code-block:: shell

   git remote add upstream https://github.com/tardis-sn/tardis.git

The ``upstream`` name is the convention used for the main TARDIS repository. The
workflow uses the HTTPS URL for upstream because it is read-only by default for
most contributors and avoids accidental writes to the main repository.

Confirm the remotes:

.. code-block:: shell

   git remote -v show

Expected output should include:

.. code-block:: text

   upstream   https://github.com/tardis-sn/tardis.git (fetch)
   upstream   https://github.com/tardis-sn/tardis.git (push)
   origin     git@github.com:your-user-name/tardis.git (fetch)
   origin     git@github.com:your-user-name/tardis.git (push)

.. _new-how-to-guide-install-tardis-in-development-mode:

Install TARDIS In Development Mode
----------------------------------

From the repository root, run:

.. code-block:: shell

   pip install -e .

TARDIS is designed to be usable from the source tree. An editable install makes
``tardis`` import from your clone, so local edits are available immediately in new
Python sessions.

After installation, confirm that Python imports from your checkout:

.. code-block:: shell

   python -c "import tardis; print(tardis.__file__)"

.. _new-tutorial-prepare-your-local-development-fork:

Prepare Your Local Development Fork
-----------------------------------

This path prepares a new contributor to work on TARDIS locally.

1. Set up a Python environment. The original developer workflow recommends
   Anaconda and refers readers to the installation guide.
2. Create a GitHub account if you do not already have one.
3. Configure your GitHub account for write access, including SSH keys.
4. Fork the TARDIS repository from the TARDIS GitHub home page.
5. Clone your fork:

   .. code-block:: shell

      git clone git@github.com:your-user-name/tardis.git
      cd tardis

6. Inspect local and remote branches:

   .. code-block:: shell

      git branch -a
      git remote -v

7. Add the main TARDIS repository as ``upstream``:

   .. code-block:: shell

      git remote add upstream https://github.com/tardis-sn/tardis.git

8. Confirm the remote configuration:

   .. code-block:: shell

      git remote -v show

   Expected remotes include ``origin``, pointing to your fork, and ``upstream``,
   pointing to ``https://github.com/tardis-sn/tardis.git``.

9. Install TARDIS in development mode:

   .. code-block:: shell

      pip install -e .

   This installs TARDIS so imports use your repository clone regardless of your
   working directory. Edits in your clone are available the next time you start a
   Python interpreter and ``import tardis``.

.. _new-how-to-guide-run-ruff:

Run Ruff
--------

TARDIS follows PEP 8 and uses Ruff for linting and formatting.

.. note::

   Pre-commit hooks are no longer used in the TARDIS developer workflow. Run
   Ruff and tests directly from the active TARDIS development environment.

Install Ruff:

.. code-block:: shell

   conda install -c conda-forge ruff

Lint code:

.. code-block:: shell

   ruff check <source_file_or_directory>

Lint and fix automatically fixable issues:

.. code-block:: shell

   ruff check <source_file_or_directory> --fix

For example, lint one parser module:

.. code-block:: shell

   ruff check tardis/io/model/parse_density_configuration.py

Or lint and fix a package area:

.. code-block:: shell

   ruff check tardis/io/model --fix

TARDIS adopts linting rules used by Astropy. Permanent rules are defined in
``pyproject.toml``. Non-permanent rules are defined in ``.ruff.toml``. Add new
rules to ``.ruff.toml``. To add a permanent rule, open a pull request against
``pyproject.toml``.

Use Ruff In VS Code
~~~~~~~~~~~~~~~~~~~

Use the official Ruff extension with the same environment you use for TARDIS:

1. Activate the TARDIS development environment in a terminal.
2. Install Ruff in that environment if it is not already installed:

   .. code-block:: shell

      conda install -c conda-forge ruff

3. In VS Code, install the Ruff extension by searching for ``Ruff`` in the
   Extensions view. You can also install it from the command line:

   .. code-block:: shell

      code --install-extension charliermarsh.ruff

4. Select the TARDIS Python interpreter with ``Python: Select Interpreter``.
   Choose the interpreter from the active TARDIS conda environment, not the base
   conda environment.
5. Add or update ``.vscode/settings.json`` in your TARDIS checkout:

   .. code-block:: json

      {
        "ruff.enable": true,
        "ruff.lint.enable": true,
        "ruff.format.enable": true,
        "editor.formatOnSave": true,
        "[python]": {
          "editor.defaultFormatter": "charliermarsh.ruff"
        },
        "editor.codeActionsOnSave": {
          "source.fixAll.ruff": "explicit",
          "source.organizeImports.ruff": "explicit"
        }
      }

With this setup, VS Code can show Ruff lint diagnostics in the editor, format
Python files on save, and apply automatic Ruff fixes when requested.

Use Ruff In Anaconda-Based Editors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Anaconda Navigator, Spyder, JupyterLab, or another editor launched from
Anaconda, make sure the editor is using the TARDIS development environment:

1. Activate the TARDIS conda environment before launching the editor, or select
   that environment in the editor's interpreter/kernel settings.
2. Confirm that Ruff is installed in that same environment:

   .. code-block:: shell

      ruff --version

3. If the editor has no Ruff integration, run the CLI commands from the
   integrated terminal:

   .. code-block:: shell

      ruff check tardis
      ruff format tardis

Avoid configuring an IDE to use Ruff from ``base`` while running tests from the
TARDIS environment, because the editor diagnostics can then disagree with the
commands used for review.
