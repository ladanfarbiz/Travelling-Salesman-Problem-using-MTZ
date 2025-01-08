Traveling Salesman Problem (TSP) Analysis

Overview

This project implements and analyzes the Traveling Salesman Problem (TSP) using the Miller-Tucker-Zemlin (MTZ) formulation. The project involves solving the TSP with exact methods using the GLPK toolkit and Python, and evaluating the runtime performance for various problem sizes.

Project Objectives
	•	Model the TSP using Mixed Integer Programming (MIP).
	•	Implement the MTZ formulation to solve the problem.
	•	Analyze the runtime performance as a function of the number of cities.
	•	Extrapolate runtime predictions for larger problem sizes.
	•	Demonstrate the limitations of exact methods and discuss heuristic alternatives.

Features
	•	Modeling: MTZ formulation with subtour elimination constraints.
	•	Runtime Analysis: Measure CPU time for varying problem sizes.
	•	Data Integration: Uses distance matrices for real-world city data.
	•	Algorithm Complexity: Analyze and model the factorial growth of runtime.
	•	Visualization: Plot runtime as a function of problem size (optional).

Technologies Used
	•	Programming Languages:
	•	Python (for MTZ implementation and runtime analysis)
	•	GNU MathProg (for modeling TSP instances)
	•	Toolkits:
	•	GLPK (GNU Linear Programming Kit) for solving MIP formulations.
	•	PuLP library for linear programming in Python.
	•	Development Tools:
	•	VirtualBox (for setting up the GLPK environment).
	•	Linux-based virtual machine (Mageia).

Project Setup

1. Install GLPK Toolkit

Follow the instructions to install the GLPK toolkit and set up a virtual machine:
	1.	Install VirtualBox.
	2.	Import the provided mip-solver.ova file into VirtualBox.
	3.	Configure a shared folder between the host and the virtual machine.
	4.	Start the virtual machine and mount the shared folder.

2. Set Up Python Environment

Install required Python libraries:

pip install pulp numpy matplotlib

3. Data File

The project uses a distance matrix (distances.dst) representing intercity distances. Ensure this file is placed in the correct directory for processing.
[ here I uploaded Cvs file]
How to Run
	1.	Solving TSP with GLPK:
	•	Load the .mod and .dat files in GLPK.
	•	Execute the following command:

glpsol -m tsp.mod -d tsp.dat


	2.	Python Implementation:
	•	Run the Python script:

python tsp_solver.py


	3.	Runtime Analysis:
	•	Modify the number of cities in the dataset.
	•	Record the runtime for each instance.
	•	Plot the results (optional).

Results

Runtime Performance

The following results summarize the CPU time required to solve TSP instances of varying sizes using the MTZ formulation:

Number of Cities	Total Distance	CPU Time (Seconds)
5	~150	0.0233
10	~550	0.2000
15	~900	9.0
20	~1200	300.0
[ these are example of results, you can adjust based on your real data and run time different on different systems]

Runtime Estimation Formula

The runtime grows factorially as the number of cities increases:

T(n) = k \times n!

Where:
	•	k = 5.51 \times 10^{-8}, derived from empirical data.

Code Structure
	•	tsp_solver.py: Python implementation of the MTZ formulation using PuLP.
	•	tsp.mod: GLPK model file defining the TSP formulation.
	•	tsp.dat: Data file containing intercity distances.
	•	City_Distance_Matrix.csv: Real-world city distance dataset.

Future Work
	•	Implement heuristic methods (e.g., genetic algorithms, nearest neighbor).
	•	Visualize the TSP solution paths for better interpretation.
	•	Optimize large instances using parallel computing.

References
	1.	Applegate et al.: Techniques for TSP solving.
	2.	Pataki: Subtour elimination in integer programming.
	3.	Wikipedia - TSP

Contact

Developed by Ladan Farbiz
Email: ladanfarbiz@yahoo.com
GitHub: https://github.com/ladanfarbiz
