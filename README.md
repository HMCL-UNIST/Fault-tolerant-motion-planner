# Fault-tolerant-motion-planner
A fault tolerant motion planner for the quadrotor subject to a complete rotor loss and rotor saturation constraints.

### Paper
Y. Nam, K. Lee, H. Shin, and C. Kwon “Fault tolerant motion planning for a quadrotor subject to complete rotor failure,” Aerospace Science and Technology, Vol. 154, Page 109529, November 2024 [(Link)](https://www.sciencedirect.com/science/article/pii/S127096382400659X)

### Directory Stucture
    Fault-Tolerant-Motion-Planner/
    │

    ├── Documents/
    
    │   ├── Full_Manuscript_With_Appendices.pdf   # Full documentation including appendices

    │   └── Appendices_Only.pdf                   # Appendices section only

    │

    └── Code/

        ├── Scene1/    # Code for scenario 1
    
        ├── Scene2/    # Code for scenario 2
    
        ├── Scene3/    # Code for scenario 3
    
        └── FTMP_Simulation_Results.m # Main file
    
   
### Dependencies
This project requires the following software:
- MATLAB 2022a
- Simulink
  
Ensure you have the correct versions installed to avoid compatibility issues.

### Usage
To run the simulations:
1. Clone the repository:
    
    - git clone https://github.com/HMCL-UNIST/Fault-Tolerant-Motion-Planner.git
    
    - cd Fault-Tolerant-Motion-Planner/Code

2. Run the Matlab code:
    - ./FTMP_Simulation_Results.m
