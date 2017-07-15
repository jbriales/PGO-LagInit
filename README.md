# PGO-LagInit
An approach for Pose Graph Optimization (PGO) initialization through the (Lagrangian) SDP relaxation.

The Semidefinite Program (SDP) relaxation for Pose Graph Optimization (or Rotation Synchronization)
may serve to recover the globally optimal solution of this problem.
However, under challenging scenarios with high rotation noise the relaxation may get suboptimal.
In those cases it is still possible to exploit the SDP solution to obtain a very good initialization.

In this repository we are releasing the *Metric Upgrade* approach,
which given a basis *V* for the primal SDP solution produces an as-feasible-as-possible candidate.
Note the Metric Upgrade approach operates on the SDP solution,
so the SDP relaxation should be solved in advance.

For details on the Metric Upgrade approach check the corresponding paper:
J. Briales and J. Gonzalez-Jimenez,
"Initialization of 3D Pose Graph Optimization using Lagrangian duality"
in Proc. IEEE Int. Conf. on Robotics and Automation (ICRA), Singapore, Malaysia, 2017.

For details on the SDP relaxation for PGO, check:
J. Briales and J. Gonzalez-Jimenez,
"Fast Global Optimality Verification in 3D SLAM"
in Proc. IEEE Int. Conf. Intell. Robots Syst. (IROS), Daejeon, South Korea, 2016.

For a very fast approach to solve the SDP relaxation, check:
J. Briales and J. Gonzalez-Jimenez,
"Cartan-Sync: Fast and Global SE(d)-Synchronization"
in IEEE Robotics and Automation Letters (RAL), 2017.
To be presented in IROS 2017, Vancouver.
Code available at https://bitbucket.org/jesusbriales/cartan-sync


## Dependencies

The current code has the following external dependencies:

- [mMath toolbox](https://github.com/jbriales/mMath) for some matrix-related operations

- [manopt toolbox](http://www.manopt.org/) for on-manifold optimization
