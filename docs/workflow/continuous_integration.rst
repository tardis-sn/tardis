**********************
Continuous Integration
**********************

We use a so-called `Continuous Integration`_ workflow with TARDIS. 
This means that each time a change is proposed (via pull-request) 
or a proposed change is merged into the main master branch a service will download the version 
and execute all the unit tests and integration tests offered with TARDIS. This helps 
us detect bugs immediately. 

The following pages explain how we setup automated debugging on TARDIS through a remote 
cloud service, called Azure, hosted by Visual Studio Team Services. This is done by 
testing, building, and securely deploying this documentation to gh-pages. As 
a developer, one should be familiar with how to smoothly run tests and record the 
documentation. 


Testing and Documentation
=========================

Testing
-------

To test TARDIS, we activate the designated environment, fetch the most recent reference data, install
TARDIS, and use the python app tester, pytest::

    $ sh ci-helpers/fetch_reference_data.sh
    $ source activate tardis
    $ python setup.py build_ext --inplace
    $ pytest tardis --tardis-refdata=$(ref.data.home) --test-run-title="TARDIS test" --cov=tardis --cov-report=xml --cov-report=html

The build_ext command is for accomdating C/C++ extentions. The --inplace command is to ensure that the build is in the source directory.

See https://docs.pytest.org/en/latest/ for more information on pytest.

Documentation
-------------

The documentation is built using Sphinx from the master branch. We use Restructured Text (reST) to encode the 
required files, which are then converted into html links, that are then pushed to the root of the 
gh-pages branch. 

See https://pythonhosted.org/an_example_pypi_project/sphinx.html and https://help.github.com/en/categories/github-pages-basics
for more information.

Setting up a secure pathway to github
=====================================

.. _install-ssh-key:

Installing a ssh-key
--------------------

When setting up the Azure pipeline for the first time, generate a ssh key locally::

    $ ssh-keygen -t rsa -b 4096 -C "your_email@example.com"

Follow along with the command prompt, and you can neglect adding a passphrase due to Azure's added security.

For more details, see `git key install`_.

Adding it as a deploy key
-------------------------

After you install the key, you must add it as a deploy key to grant access to push to your repository.

To do this, go to GitHub_
    - Go to your profile
    - Go to your desired repository
    - Open the Settings
    - Click on Deploy keys and copy the contents of your public key (The default name is id_rsa.pub).
                
For explicit details, see `git deploy key`_.

Adding your key locally and copying the known host name
-------------------------------------------------------

If you wish to deploy the documentation locally to gh-pages, you must add the generated key on your computer::

    $ eval "$(ssh-agent -s)"
    $ ssh-add ~/.ssh/id_rsa (Or whatever you called your key)
    $ ssh-agent -k

Copy the saved known host, as we will need it for installing the key on Azure.
It should look something like (should look something like [1]As3..=ssh-rsa ..) and will be the last line added to::

     ~/.ssh/known_hosts 


Setting up Azure services 
=========================

The first step is to visit `Azure Devops`_ and create an account.

Adding the key to Azure's secure files
--------------------------------------

Upload the secure key to Azure's library to download it in the script.

In the library tab, select the **Secure files** tab, upload the file, and authorize it for use in all pipelines.
As seen here:
  .. image:: images/secure_file.png

See also: `Azure secure files`_   

One must use the DownloadSecureFile@1 `Azure task`_ to download the file onto the virtual machine when the pipeline runs.
After which, one must use the InstallSSHKey@0 `Azure task`_ to add the ssh key.

Setting up the YAML file to deploy
----------------------------------

YAML is short for "YAML Ain't Markup Language", as it is intented to be a simple way to write script
that is standard for all programing languages. It is the file that communicates directly with the
pipeline. 

In the script, specify which branches you want to trigger for continuous deployment 
and/or applicable for pull requests. Otherwise, it triggers all of them::

    trigger:
    - branch_name

    pr:
    - branch_name

It follows the following hierarchy:

Pipeline hierarchy
^^^^^^^^^^^^^^^^^^

.. graphviz::

  digraph {
    a -> b -> c -> d -> e
    c -> f
    a [label="Stages",shape=circle,fillcolor="white",style="filled"];
    b [label="Stage",shape=circle,fillcolor="white",style="filled"];
    c [label="Jobs", shape="circle", fillcolor="white", style="filled"]
    d [label="Job", shape="circle", fillcolor="white", style="filled"]
    e [label="Steps", shape="circle", fillcolor="white", style="filled"]
    f [label="Task", shape="circle", fillcolor="white", style="filled"]
  }

A job is referred to as the most basic building block that the pipeline runs, using a single agent, which 
is composed of steps of script.
A task is a predefined script with a definitive purpose, such as downloading the secure file or installing a ssh-key.

To download and add the ssh-key, prepare the scripts as::

      - task: DownloadSecureFile@1
        inputs: 
          secureFile: 'id_azure_rsa'

Secure files stored in the Azure server are encryped and again decrypted by the Azure task that uses the file.

Download a secure file to a temporary location in the virtual machine::

      - task: InstallSSHKey@0
        inputs:
          hostName: $(gh_host)
          sshPublicKey: $(public_key)
          #sshPassphrase: # Optional - leave empty if it was left empty while generating the key.
          sshKeySecureFile: 'id_azure_rsa'


hostName is the line that was copied in `Adding your key locally and copying the known host name`_.

sshPublicKey should be a string value of what is inside your .pub file (i.e: rsa-key Axddd... username@server).

sshKeySecureFile is the downloaded secure file you generated, you can reference directly as shown.

For more details, see `Azure ssh-task`_

To define variables in the script, one can do so using key-value pairs. For example::

    variables:
      system.debug: 'true'

To define secret variables, or variables outside the script, one must navigate to variables after
selecting the three dots on the top right while editing that pipeline, as seen here:  

  .. image:: images/variables.png

After defining the variable, one can encrypt it using this lock symbol:

  .. image:: images/lock.png

Variables are referenced as $(variable_name), as seen in the InstallSSHKey@0 task in the hostName and sshPublicKey inputs.

Azure provides a list of agent hosts that can run the pipeline on a virtual machine. In our pipelines, we
use the vm_Images: Ubuntu 16.04 and macOs-10.13.

In a job, you can list a single vm as::

      pool:
        vmImage: "Ubuntu 16.04"

If you are using a self hosted agent (see `Installing and running a self hosted agent` for more details)::

      pool:
        name: "agent_pool_name"

Or if you prefer to use multiple virtual machines and specity the maximum that can run at the same time::

      strategy:
        matrix: 
          linux:
            vmImage: "Ubuntu 16.04"
            conda: '/usr/share/miniconda'
          mac:
            vm_Image: 'macOS-10.13'
            miniconda.url: 'http://repo.continuum.io/miniconda/Miniconda2-latest-mac-x86_64.sh'
        maxParallel: 4



Additional references
--------------------

https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html#inline-markup-and-special-characters-e-g-bold-italic-verbatim

https://docs.microsoft.com/en-us/azure/devops/pipelines/?view=azure-devops

https://yaml.org/

https://docs.microsoft.com/en-us/azure/devops/pipelines/yaml-schema?view=azure-devops&tabs=schema


.. include:: git_links.inc
.. include:: azure_links.inc