GRINS
=======

General Reacting Incompressible Navier-Stokes (GRINS) was initiated                                                                                                                                                      
to house common modeling work centered around using the incompressible
and variable-density (low-Mach) Navier-Stokes equations
utilizing the [libMesh](https://github.com/libMesh/libmesh.git) finite
element library, including both MPI and MPI+threads parallelism,
as provided by [libMesh](https://github.com/libMesh/libmesh.git). 
GRINS has now become a tool for rapidly developing
formulations and algorithms for the solution of complex multiphysics
applications. 
GRINS originally lived within
the [PECOS](http://pecos.ices.utexas.edu) center at the Institute for Computational
Engineering and Sciences ([ICES](https://www.ices.utexas.edu))
at [The University of Texas at Austin](https://www.utexas.edu).

We encourage pull requests for new features, bug fixes, etc. For questions regarding development,
we have a grins-devel Google group setup. For user related questions, please use the grins-users
group.

Dependencies
============

Requirements
------------

In addition to a modern C++ compiler,
GRINS requires an up-to-date installation of the [libMesh](https://github.com/libMesh/libmesh.git)
finite element library (currently version 0.9.2.2 or later) as well as the Boost C++ library.

Optional Packages
-----------------

To enable the reacting low Mach Navier-Stokes physics class, GRINS must be compiled with
an external chemistry library. While [Cantera](http://code.google.com/p/cantera/) is
partially supported, [Antioch](https://github.com/libantioch/antioch) is fully
supported.

Building GRINS 
================

GRINS uses an Autotools build system, so typical GNU build commands are used.

1. ./bootstrap (generate configure script)
2. ./configure --prefix=/path/to/install --with-libmesh=/path/to/libMesh --with-boost=/path/to/boost (for more options, do ./configure --help)
3. make (note parallel builds are supported)
4. make check (note parallel-tests are supported)
5. make install

METHOD
------

By default, GRINS leverages the METHOD environment variable
(described [here](https://github.com/libMesh/libmesh/blob/master/README.md)) in order to
retrieve the CXXFLAGS variable from the [libMesh](https://github.com/libMesh/libmesh.git)
installation (if METHOD is not present, the default is "opt"). The user can define
their own CXXFLAGS variable by passing 
<pre><code>
--disable-libmesh-flags CXXFLAGS="your flags here"
</code>
</pre>
to configure.
