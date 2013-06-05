
=== how to run this code ===

  To run iPic3D on the Xeon Phi processor (knc2-mic0):

    ssh miclogin
    ssh knc2
    cd iPic3D
    make -f makefile.mic clean
    make -f makefile.mic -j3
    make -f makefile.mic run

  To run iPic3D on the Xeon processor (knc2):

    ssh miclogin
    ssh knc2
    cd iPic3D
    make -f makefile.knc clean
    make -f makefile.knc -j3
    make -f makefile.knc run

  To modify the code in /home/eajohnson/iPic3D:
  * Edit the "run" target of the makefile to change the number
    of MPI processes (e.g. "-n 16") and the number of OMP thread per process.
    If you change the number of MPI processes, then you must:
    + Modify "processtopology/VCtopology3D.h" (search for "XLEN =")
      so that XLEN times YLEN times ZLEN equals number of MPI processes.
      (For the 2D test problem I keep ZLEN=1.")
    + Recompile (make clean and make).
    + If necessary, edit "inputfiles/Random.inp",
      making sure that XLEN divides nxc, YLEN divides nyc,
      and ZLEN divides nzc.
