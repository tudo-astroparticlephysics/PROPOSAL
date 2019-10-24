Installation
------------

PRopagate with Optimal Percision and Optimal Speed for All Leptons works on the
basis of C++ to maximize the speed of propagation.
To install PROPOSAL with pyBindings you have to install the dependencies first.
An example for Ubuntu 16.04 is given below.

.. code-block:: shell

   apt install cmake doxygen liblog4cplus-dev

Package names may vary depending on OS.
Then you can install PROPOSAL at your preferred location.
Create a build and source directory

.. code-block:: shell

   mkdir -p PROPOSAL/{build,src}

and clone the `repo <https://github.com/tudo-astroparticlephysics/PROPOSAL>`_ from github

.. code-block:: shell

   git clone https://github.com/tudo-astroparticlephysics/PROPOSAL.git PROPOSAL/src

Move to the build directory and generate the Makefile with cmake.

.. code-block:: shell

   cd PROPOSAL/build
   cmake ../src \
      -DADD_PYTHON=ON \
      -DADD_PERFORMANCE_TEST=OFF \
      -DADD_ROOT=ON \
      -DADD_TESTS=OFF \
      -DCMAKE_INSTALL_PREFIX=/custom/prefix \
      -DPYTHON_LIBRARY=$(python-config --prefix)/lib/libpython2.7.dylib \
      -DPYTHON_INCLUDE_DIR=$(python-config --prefix)/include/python2.7

Compile the project and install it

.. code-block:: shell

   make -j && make install

Don't forget to add your install directory to your .bashrc

.. code-block:: shell

   export PYTHONPATH=${PYTHONPATH}:/custom/prefix/lib
   export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:custom/prefix/lib
