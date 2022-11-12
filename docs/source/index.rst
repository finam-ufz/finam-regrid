:html_theme.sidebar_secondary.remove: true

============
FINAM Regrid
============

Regridding adapter for the `FINAM <https://finam.pages.ufz.de/>`_ model coupling framework,
using `ESMPy <https://earthsystemmodeling.org/esmpy/>`_.

Quickstart
----------

Installation:

The easiest way to install ESMPy is to use `Conda <https://docs.conda.io>`_.

Create a Conda environment with ESMPy installed, and activate it:

.. code-block:: Shell

    $ conda create -n esmpy -c conda-forge esmpy
    $ conda activate esmpy

For details, see the ESMPy `documentation <https://earthsystemmodeling.org/esmpy_doc/release/latest/html/>`_.

Now, install :mod:`finam_regrid` and other required packages using `pip <https://pip.pypa.io>`_:

.. code-block:: Shell

    $ pip install git+https://git.ufz.de/FINAM/finam-regrid.git

Usage
-----

The main class of this package is the adapter :class:`.Regrid`.
For other package members, see the :doc:`api/index`.

FINAM-regrid supports all of ESMPy's  :class:`regridding methods <.RegridMethod>`.
For more options, see the ESMPy class `Regrid <https://earthsystemmodeling.org/esmpy_doc/release/latest/html/regrid.html>`_

See also the `example scripts <https://git.ufz.de/FINAM/finam-regrid/-/tree/main/examples>`_
in the GitLab repository for fully functional usage examples.

API References
--------------

Information about the API of FINAM-regrid.

.. toctree::
    :hidden:
    :maxdepth: 1

    self

.. toctree::
    :maxdepth: 1

    api/index
