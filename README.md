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
we have a [grins-devel](https://groups.google.com/forum/#!forum/grins-devel) Google group setup. For user related questions, please use the [grins-users](https://groups.google.com/forum/#!forum/grins-users)
group.

Dependencies
============

Requirements
------------

In addition to a C++11 compiler, GRINS requires an up-to-date installation of the [libMesh](https://github.com/libMesh/libmesh.git) finite element library. 

libMesh
-------
GRINS development both drives and is driven by libMesh development. Thus, the required minimum master hash of libMesh may change in GRINS master. The current required libMesh master hash is 2d9b8bc, as of GRINS [PR #520](https://github.com/grinsfem/grins/pull/528).
GRINS release 0.5.0 can use libMesh versions as old as 0.9.4. Subsequent to
the 0.5.0 release requires at least libMesh 1.0.0.


Optional Packages
-----------------

To enable the reacting low Mach Navier-Stokes physics class, GRINS must be compiled with
an external chemistry library. Both [Cantera](http://code.google.com/p/cantera/) and
[Antioch](https://github.com/libantioch/antioch) are fully supported.

The current required minimum hash for using Antioch is libantioch/antioch@e17822d (libantioch/antioch#265).


Building GRINS 
================

GRINS uses an Autotools build system, so typical GNU build commands are used. We support, and encourage, out-of-source builds (so-called VPATH builds). However, in-source builds are supported.

1. cd grins-clone
2. ./bootstrap (generate configure script)
3. cd ../ & mkdir build & cd build
4. ../grins-clone/configure --prefix=/path/to/install --with-libmesh=/path/to/libMesh (for more options, do ../grins-clone/configure --help)
5. make (note parallel builds are supported)
6. make check (note parallel-tests are supported)
7. make install

LD_LIBRARY_PATH
---------------

If you've compiled libMesh with PETSc or other external libraries and have compiled GRINS with Antioch, Cantera, or other external libraries, you will need to add them to your LD_LIBRARY_PATH as we do not use -rpath when linking to the libraries.

METHOD
------

By default, GRINS leverages the METHOD environment variable
(described [here](https://github.com/libMesh/libmesh/blob/master/README.md)) in order to
retrieve the CXXFLAGS variable from the [libMesh](https://github.com/libMesh/libmesh.git)
installation (if METHOD is not present, the default is "opt"). Note that unlike libMesh,
GRINS currently only supports building one METHOD at a time. Hence, we use `METHOD` and
not `METHODS`. For example
<pre><code>
./configure METHOD=devel
</code>
</pre>
is valid.

The user can define
their own CXXFLAGS variable by passing 
<pre><code>
--disable-libmesh-flags CXXFLAGS="your flags here"
</code>
</pre>
to configure.

Examples
========

Upon running `make install`, there are several examples in the `/path/to/install/examples` directory. Each example can be run with the local `run.sh` script. You may set the environment variable `GRINS_RUN` to run with more than one processor, e.g. `export GRINS_RUN="mpiexec -np 4"`. Additionally, you can set the environment variable `GRINS_SOLVER_OPTIONS` to pass solver options, e.g. to use MUMPS through PETSc (if you built libMesh with PETSc and built PETSc with MUMPS), `export GRINS_SOLVER_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps"`.
