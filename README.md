# vtkNonRigidICP
Non-Rigid ICP based off optimal step non rigid ICP paper from CVPR 2007. Implemented as a VTK Filter using Eigen for sparse matrix operations. Build with provided CMakeLists.txt.

Implements the data term and regularization term from the paper. The data term drives points to points whilst the regularization terms smoothens the transform among its neighbors.

Filter can cache the results to pick up where it left off faster. It can also output the displacement vectors as scalars.

Future improvements: Better correspodence dropping. Incoporate normals in kd-trees nearest neighbor search.

![Alt text](https://andaharoo.files.wordpress.com/2019/09/screenshot.png)