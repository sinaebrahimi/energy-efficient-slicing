# energy-efficient-slicing
MATLAB code for the simulation of my paper entitled "Joint Resource and Admission Management for Slice-enabled Networks"

S. Ebrahimi, A. Zakeri, B. Akbari, N. Mokari, "Joint Resource and Admission Management for Slice-enabled Networks",  IEEE/IFIP Network Operations and Management Symposium (NOMS 2020), Budapest, Hungary, 2020.

IEEE Xplore:
https://ieeexplore.ieee.org/document/9110322

ResearchGate:
https://www.researchgate.net/publication/337673407_Joint_Resource_and_Admission_Management_for_Slice-enabled_Networks

The code includes two joint and disjoint methods which are compared in the paper. The execution/reading the code should be commenced from energy_efficient_slicing-*.m files.

The two optimization problems include an admission control before solving the optimization problem. We used the MOSEK solver and CVX package to solve problems (8), (10), (11), (12), (13), and (14).
Moreover, all simulation steps (including initialization, admission control mechanisms, and solving ILP problems with MOSEK toolbox) have been implemented in MATLAB software which is widely used to solve resource allocation problems.

Network architecture of the CSP in our system model:
![System Model (Architecture)](https://github.com/sinaebrahimi/energy-efficient-slicing/blob/master/figures/1-Network%20architecture%20of%20the%20CSP%20in%20our%20system%20model-transparent.png?raw=true "System Model (Architecture)")

Simplified flowchart of JRA (Joint Resource Allocation method):
![JRA Flowchart](https://github.com/sinaebrahimi/energy-efficient-slicing/blob/master/figures/2-Simplified%20flowchart%20of%20JRA%20(Joint%20Resource%20Allocation).png?raw=true)

Simplified flowchart of DRA (Disjoint Resource Allocation method):
![DRA Flowchart](https://github.com/sinaebrahimi/energy-efficient-slicing/blob/master/figures/(Not)-Simplified%20flowchart%20of%20DRA%20(Disjoint%20Resource%20Allocation).png?raw=true)

Comparing the total cost between DRA and JRA methods:
![Total Cost Comparison](https://github.com/sinaebrahimi/energy-efficient-slicing/blob/master/figures/3a-Overall%20Cost%20(C_total).png?raw=true)

To get more useful figures including flowcharts and graphs refer to the 'figures' folder. Some of them are not included in the paper.
