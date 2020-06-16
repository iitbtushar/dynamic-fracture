# dynamic-fracture
Reference: Evaluation of variational phase-field models for dynamic brittle fracture
-- by T.K. Mandal, V.P. Nguyen and J.Y. Wu, Engineering Fracture Mechanics (June, 2020)

The repository contains the main input and output data for the above-mentioned paper. This includes following directories correspond to different dynamic fracture examples using PFM as described in the paper. FE mesh have been generated using GMSH and simulations are carried out using in-house C++ based FE code feFRAC based on the open source jive library at https://software.dynaflow.com/jive/

1. "branching": Dynamic crack branching for a Rectangular strip subjected to tensile impulse traction (Sec. 3.1);
2. "branching2": Rectangular strip subjected to tensile impulse traction on the pre-crack faces (Sec. 3.2);
3. "fragmentation": Dynamic fragmentation of a thick cylinder (Sec. 3.3);

4. "arrest": Dynamic crack arrest due to holes (Sec. 4.1);
5. "velocity-toughening": Crack branching and velocity-toughening mechanism (Sec. 4.2);
6. "multiple-branching": Dynamic crack arrest due to holesMultiple branching (Sec. 4.3);