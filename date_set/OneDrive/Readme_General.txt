Authors: 
	Francesco Niosi, 	email: francesco.niosi@polito.it 
	Oronzo Dell'Edera, 	email: oronzo.delledera@polito.it
	Sebastian Schreier,	email: s.schreier@tudelft.nl	
Credit: Technology University of Delft
Date: 19/10/23
Last Edit: 20/10/23

###########################################################################################################################################################################################

In this folder, you will find all the experimental tests conducted during the experimental campaign described in the paper "Experimental Investigation and Dataset Release of a Taut-Leg Mooring System for a Semi-Submersible Floating Offshore Wind Turbine."

###########################################################################################################################################################################################

-The "ALS" folder contains tests in the Accidental Limit State where the turbine and/or mooring system is damaged.

-The "EmptyTank" folder contains all the wave tests categorized by wave type (Irregular, Regular, System ID). 

-The "FreeDecays" folder contains all the FreeDecay tests with and without mooring.

-The "LCS_DesignPretension" folder contains tests performed with the spring configuration, where the lines are tensioned at the design pretension level.These tests are also categorized by wave type. 

-The "OffDesignPretension" folder contains tests with pretension levels different from the design, and it is divided into LCS and LCN, with LCS referring to the spring configuration and LCN to nylon.

###########################################################################################################################################################################################

All the data are acquired at Fs=100Hz and filtered at 10Hz. The translational motions (Surge, Sway, Heave) are referenced to the center of gravity (CoG) of the floating structure.
Detailed information about the organization of each file can be found in the readme inside each folder.

###########################################################################################################################################################################################

Description of monitored channels:
-Time			=	Time vector
-Surge			=	Surge time history
-Sway			=	Sway time history
-Roll			=	Roll time history
-Pitch			=	Pitch time history
-Yaw			=	Yaw time history
-MooringLine1		=	Loadcell acquired signal time history of mooring line 1
-MooringLine2		=	Loadcell acquired signal time history of mooring line 2
-MooringLine3		=	Loadcell acquired signal time history of mooring line 3
-WindForce		=	Loadcell acquired signal time history of wind winch
-WaveProbe1		=	Wave Gauge 1 time history
-WaveProbe2		=	Wave Gauge 2 time history
-WaveProbe3		=	Wave Gauge 3 time history
-WaveProbe4		=	Wave Gauge 4 time history
-WaveProbe5		=	Wave Gauge 5 time history
-WMExtPosition		=	Wave Maker position time history with respect to its rest position. AGGIUNGERE POSITIVE NEGATIVE
-WMSteering		=	Wave Maker input in voltage
-MarkerError		=	Optitrack error time history
-Trigger		=	Trigger signal

###########################################################################################################################################################################################