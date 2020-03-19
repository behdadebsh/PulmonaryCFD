# PulmonaryCFD
This repository includes the simulations and setups for 3D CFD simulations

There are three main directories included here:
1) SingleBif ---> This is a case for a single bifurcation which the mesh for it was created using ICEMCFD (Ansys). This was coverted to OpenFoam (OF) mesh within OF and a steady-state simulation was performed with constant boundary conditions.
2) Bifurcations ----> This has a mesh setup using a surface stl for a mesh having the 3rd generation bifurcations (4 outlet branches). The mesh was created using snappyHexMesh and the results are within the folder.
3) BIF_CFMesh ----> Using a perpendicular capped stl file similar to Bifurcations case. Although the mesh has been created using CfMesh for simulations. This case has the best mesh quality among all cases mentioned. This case has the option of showing residuals (for following the convergence while simulation is running - the functions is included in controlDict). Also, it has the parabolic velocity profile function in U (velocity) initial boundary condition but currently is commented out since the plots in the case folder are relevant to uniform inlet BC.

