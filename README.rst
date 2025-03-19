==================================
Schrodinger Maestro Scipion plugin
==================================

**Documentation under development, sorry for the inconvenience**

This is a **Scipion** plugin that offers different **Schrodinger tools**.
These tools will make it possible to carry out different functions for Virtual Drug Screening (VDS)

Therefore, this plugin allows to use programs from the Schrodinger software suite
within the Scipion framework.

Full documentation to this plugin can be found in the `official documentation page <https://scipion-chem.github.io/docs/plugins/schrodinger/index.html>`_.

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**   (currently not available)

.. code-block:: 

      scipion installp -p scipion-chem-schrodinger
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem-schrodinger.git

2. **Switch to the desired branch** (master or devel):

Scipion-chem-schrodinger is constantly under development and including new features.
If you want a relatively older an more stable version, use master branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-schrodinger
            git checkout devel

3. **Install**:

.. code-block::

            scipion3 installp -p path_to_scipion-chem-schrodinger --devel

- **Binary files** 

The **Schrodinger software must be installed separately**, as it requires a license.
You must either install it in the EM_ROOT (typically: SCIPION_HOME/software/em/Schrodinger202*-*)
or define the Schrodinger path in scipion.conf as SCHRODINGER_HOME

- **Tests**

To check the installation, simply run the following Scipion test:

scipion3 tests schrodinger.tests.test_sitemap.TestSitemap

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/schrodinger_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/schrodinger_prod.svg
