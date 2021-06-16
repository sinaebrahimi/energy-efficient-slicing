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

![System Model (Architecture)](https://github.com/[username]/[reponame]/blob/[branch]/figures/1-Network%20architecture%20of%20the%20CSP%20in%20our%20system%20model.png?raw=true)
