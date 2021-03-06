# Vibration analysis of a simplified model of an industrial washing machine with Matlab and Simulink.


Excitations by an out of balance force are of great importance for many technical applications. 
For example, the laundry in a washing machine can cause an unbalance, 
which in the worst case, if the motor is too small, can stop the whole machine. 
In this university course-work, we take a closer look at the unbalance excitation that occurs during the centrifuge of an industrial washing machine and try to eliminate the unwanted vibrations with the help of a compensation damper. In doing so we have taken a simplified model of a real washing machine.
From the simulations performed, we can see how well an undesirable vibration can be  can be compensated for by means of an absorber and the influence of the 
parameter selection on the overall result.

The original task definiton is: 

An industrial washing machine (Figure 1) is located in a hall and is supported on the the floor with a spring-damper system. Due to the unbalance unacceptably high vibration amplitudes occur in nominal operation (centrifuge). To solve this problem, the machine must be stabilized by means of a vibration damper. The following points are considered:
- calculation of vibration amplitudes in nominal operation without absorber.
- calculation of the required spring stiffness of the absorber.
- calculation of vibration amplitudes in nominal operation with absorber. 
- presentation of the amplitude frequency response for the damped system. 
- evaluation of the influence of the absorber stiffness on the results.

This repo contains the code used to solve the points mentioned above. It consists of the following files: 
- Unwuchterregung_Anlauf.m - Matlab script that estimates the excitation during the run-up phase of the washing machine.
- Unwuchterregung_mit_Tilger.m - Matlab script that estimates the excitation during a centrifuge stage with a vibration damper
- Unwuchterregung_ohne_ST.m - Matlab script that estimates the excitation during a centrifuge stage with a vibration damper
- Unwucht4_Frequenzgang.slx - Simulink model that simulates the frequency response of the system

## License
MIT
