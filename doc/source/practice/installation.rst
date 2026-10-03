First Simulation
----------------

This section is a complete tutorial to download and run your first **HarP** simulation.

HarP relies on a Unix-like environment (such as Linux or macOS) to run its Makefile. The Makefile creates and sets up all the libraries required to run your first simulation in a specific Python environment.

.. note::
    **For users who work with their own python environment**, you can download the repository directly from the `HarP GitHub <https://github.com/Syttam/HarP>`_ page and go directly to the **First Run** section. 
    We do not recommend this approach since problems can arise from the user libraries versions. Make sure to have *numpy, scipy, h5py* libraries up to date for the run, as well as *jupyter notebook* for the example of post-process analysis.

----

Windows Users
^^^^^^^^^^^^^

If you are on Windows, the officially supported and recommended workflow is to use the **Windows Subsystem for Linux (WSL2)**.

1. Open PowerShell as Administrator and install WSL with Ubuntu:

   .. code-block:: powershell

       wsl --install

2. Once your computer restarts and your Ubuntu terminal is ready, ensure the basic build tools and Python are installed:

   .. code-block:: bash

       sudo apt update
       sudo apt install make python3 python3-venv python3-pip

3. **Important for performance:** Ensure you clone the HarP repository directly into your Linux file system (e.g., ``~/src/``) rather than your mounted Windows drive (``/mnt/c/``). File operations are significantly faster inside the native Linux filesystem.

From this point on, you can follow the standard Linux download and installation instructions using your WSL terminal. Modern IDEs like VS Code and PyCharm feature dedicated WSL extensions that allow you to edit the code in Windows while the solver runs seamlessly in Linux.


----

Linux or Mac Users
^^^^^^^^^^^^^^^^^^

Clone the latest version of HarP from GitHub:

.. code-block:: bash

    cd path/where/you/want/to/have/it/
    git clone https://github.com/Syttam/HarP.git


HarP uses a Makefile to simplify environment setup and dependency management across different operating systems.

.. warning::
    The Makefile requires Python 3.13+ (CI tests 3.13 and 3.14).

To create an isolated Python virtual environment (``.venv``) and install the core dependencies required to run your first simulation, navigate to the repository root and run:

.. code-block:: bash

    make init

If your plan is to use the post-processing Jupyter notebook provided, you should instead install the development dependencies. 
This command will set up the environment with all necessary analysis and testing packages:

.. code-block:: bash

    make dev

Once completed, remember to activate the virtual environment before running any scripts manually:

.. code-block:: bash

    source .venv/bin/activate

You are now in position to make your first simulation :ref:`first-run`.

----

.. _first-run:

First Run
^^^^^^^^^

Before running your first simulation, it is recommended to create a dedicated working directory to store all input and output files.
This keeps your simulation files separate from the HarP source code, which generally does not need to be modified.


.. code-block:: bash

    mkdir first_simulation
    cd first_simulation


An example input file, ``input.py``, is provided in the repository. It describes a simple physical configuration in which a :math:`810` nm laser with :math:`2.5 \times 10^{14}` W/cm\ :sup:`2` peak vacuum intensity propagates through a uniform Argon gas cell.
The simulation computes the fundamental field together with harmonic orders ranging from 3 to 59, as set by the ``order`` array in ``input.py``.

Copy the example input file, ``input.py``, into your working directory:

.. code-block:: bash

    cp path/to/HarP/input.py .

By default, ``input.py`` writes its diagnostics next to itself (via ``diag_path = os.path.dirname(os.path.abspath(__file__)) + os.sep``), so output lands in your working directory rather than inside the :program:`HarP` repository - override ``diag_path`` in your own copy if you want it somewhere else.

.. warning::

    The path that points to the tables responsible for the microscopic description, ``tables_path`` (line 80 of ``input.py``), is resolved relative to the directory :program:`HarP` is *run from*. If you run the simulation from somewhere other than the :program:`HarP` repository itself, replace::

        tables_path = 'tables/'

    with the absolute path::

        tables_path = 'path/to/HarP/tables/'


The diagnostics are defined at the end of ``input.py``. By default, the example enables four diagnostics:
    * the harmonic energy as a function of the propagation distance,
    * the electric field envelope :math:`\mathcal{E}_q` in :math:`(r,z,t)` space,
    * the atomic density in :math:`(r,z,t)` space,
    * the electronic density in :math:`(r,z,t)` space,

Additional diagnostics, such as gas density profiles or far-field distributions, can be enabled before running the simulation.

You can then start the simulation with:

.. code-block:: bash

    python path/to/HarP/HarP.py input.py

.. warning::

    Make sure that you have enough memory to store the different diagnostics (:math:`\sim 300` MB of data with the provided input file.)


As an indication, using :math:`N_q = 30`, :math:`N_r = 512`, :math:`N_z = 256` and :math:`N_t = 256`, the simulation completed in approximately :math:`130` s (:math:`\sim 2` min) on a MacBook Air M5 with 16 GB of memory. 
The time scaling is linear with :math:`N_t,N_z,N_r` and :math:`N_q`.

Once the simulation is complete the different diagnostics are stored in **HDF5** files, and you can continue with the :doc:`Post process analysis <../practice/post_process>`.

.. note::

    To run the simulation, you can also copy the launcher ``HarP.py`` in your working directory and modify the 5th line such that it refers to the absolute path of the *src* folder:
    
    .. code-block:: python
        
        sys.path.append(os.path.join(os.path.dirname(__file__), 'path/to/HarP/src'))

    Then you can directly start the simulation in your working directory with:

    .. code-block:: bash

        python HarP.py input.py

----

Additional features
^^^^^^^^^^^^^^^^^^^

To access the documentation from your own laptop and not from the GitHub web pages, you can compile and open it with:

.. code-block:: bash

    make docs
    open doc/build/html/index.html


You can also verify that the installation was successful (after the ``make dev`` command) and that the code performs as expected on your hardware architecture with:

.. code-block:: bash

    make test

The test suite mainly covers computational performance so far but should also cover correctness. If you wish to run only the performance benchmarks to evaluate the solver's time complexity and memory scaling on your machine, you can instead use:

.. code-block:: bash

    make benchmark

**A Note on Continuous Integration (CI)**

You might see references to "Continuous Integration" or "CI" in the repository or on GitHub. In software development, CI is an automated quality-control system. Every time a new commit is pushed to GitHub, remote servers automatically launch fresh Linux and macOS environments, install HarP, and run the ``make test`` suite. This guarantees that new physics models or code adjustments do not accidentally break the solver for other users or on different operating systems.
