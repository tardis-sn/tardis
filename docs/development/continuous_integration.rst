**********************
Continuous Integration
**********************

We use a so-called `continuous integration`_ workflow with TARDIS.
This means that each time a change is proposed (via pull request)
or a change is merged into the *master* branch, a service will
clone the repository, checkout to the current commit and execute
all the TARDIS tests. This helps us to detect bugs immediately.


Azure Pipelines
===============

Currently, we use the `Azure DevOps`_ service to run most of our
pipelines. The following section explains briefly the different
components of a pipeline.


YAML files
----------

A pipeline is essentially a YAML configuration file with different
sections such as variables, jobs and steps. Unlike other services
such as GitHub Actions, pipelines on Azure must be created through
the web UI for the first time. Then, making changes to an existing
pipeline is as easy as making a pull request.


Triggers
--------

First thing to do is telling the pipeline when it should run.
*trigger* (also known as the CI trigger) sets up the pipeline to
run every time changes are merged to the *master* branch.
::
  trigger: 
    - master

If some trigger is not specified then the default configuration
is assumed.
::
  trigger:
    branches:
      include:
      - '*'

  pr:
    branches:
      include:
      - '*'

This means the pipeline will start running every time changes are 
merged  to any branch of the repository, or someone pushes new
commits to a pull request.

If you want to run a pipeline only manually set both triggers to 
*none*.
::
  trigger: none

  pr: none

Notice that you can test changes in a pipeline by activating the PR
trigger on a new pull request, even if that trigger is disabled on
the YAML file present in the *master* branch.

There are more useful triggers such as the *cron* trigger, see the 
`Azure documentation section on triggers`_ for more information.

.. warning:: Triggers also can be set on the Azure's web interface 
          too, but this action is discouraged, since it overrides
          any trigger specified in the YAML file and could lead to
          confusing sitations.

Variables
---------

Jobs
----

A job is referred to as the most basic building block that the pipeline runs, using a single agent, 
which is composed of steps of script.

A task is a predefined script with a definitive purpose, such as downloading the secure file or installing a ssh-key.

To download and add the ssh-key, prepare the scripts as::

      - task: DownloadSecureFile@1
        inputs: 
          secureFile: 'id_azure_rsa'

Secure files stored in the Azure server are encrypted and again decrypted by the Azure task that uses the file.

Download a secure file to a temporary location in the virtual machine::

      - task: InstallSSHKey@0
        inputs:
          knownHostsEntry: $(gh_host)
          sshPublicKey: $(public_key)
          #sshPassphrase: # Optional - leave empty if it was left empty while generating the key.
          sshKeySecureFile: 'id_azure_rsa'


In a job, you can list a single vm as::

      pool:
        vmImage: "Ubuntu 16.04"

Or if you prefer to use multiple virtual machines and specify the maximum that can run at the same time, in 
addition to specifying variables as key-value pairs such as conda and miniconda.url below. ::

      strategy:
        matrix: 
          linux:
            vmImage: "Ubuntu 16.04"
            conda: '/usr/share/miniconda'
          mac:
            vm_Image: 'macOS-10.14'
            miniconda.url: 'http://repo.continuum.io/miniconda/Miniconda2-latest-mac-x86_64.sh'
        maxParallel: 2
      pool:
        vmImage: $(imageName)

This trick is also convenient for specifying different variable builds for the same vmImage. As one can keep the vm_Image 
constant and change the key value pair for each job in the matrix.

.. include:: git_links.inc
.. include:: azure_links.inc
