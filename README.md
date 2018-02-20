# Developing a Hardware-in-the-Loop Testbed for Smart Building Energy Management Systems
Proposal Code:
-
8109

Main Supervisor (email)	
-
Peng Chih-Hsien Jimmy (elepcj{@}nus.edu.sg )

Title:
-
Developing a Hardware-in-the-Loop Testbed for Smart Building Energy Management Systems

Synopsis:	
-
In a liberalised energy market, demand response (DR) is an effective mechanism for optimising the dispatch of generations and achieve systemwide benefits. DR aggregators will play a key role in managing the DR of a collective of customers. This leads to the concept of building energy management system (BEMS), which is a software-cum-hardware tool that controls the various components of the residential distribution system, such as the loads, photovoltaic generation, battery storage, etc. While various approaches exist for designing the software component of the BEMS, the hardware implementation remains a challenge. The control circuits for each of the controllable loads and generation sources need to effectively communicate with a central BEMS controller, which is essentially software running on a dedicated microcomputer or a PC. 

Objective:
-
Develop a BEMS hardware component that will be interfaced to another PC running a power system simulation tool. The power system measurements will be input to the BEMS controller, and the outputs of the BEMS controller will be fed back to the virtual power system. 

Nature of Work:
-
Modelling, Simulation, Embedded system 

Project deliverables:
-
Development of BEMS controller on an appropriate computation platform such as Arduino, DSP, or others. The BEMS controller should communicate with a PC running power system simulation software, such as MATLAB dSPACE. The designed BEMS controller should be flexible so as to implement and test any smart grid functions as required.

----------------------------------------------------Repository Description----------------------------------------------------
-
The MATLAB functions NL_LP, IL_LP, IL_EM_LP, TCL_LP will be called to solve the optimization problem of minimizing cost of electricity for non-interruptible, interruptible without electrical machinery, interruptible with electrical machinery and thermostatically controlled loads respectively. These functions will optimize the objective function by using the 'linprog' function.

[updated 20/02/18]
