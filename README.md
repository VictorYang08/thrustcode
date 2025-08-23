# thrustcode
This repository houses the code behind my paper on optimizing Northrop Grumman's GEM63XL. It includes the code for ideal and applied performance. 

# IMPORTANT NOTES
The ideal performance code is housed under ```thrustvtime_v2.py``` and the applied performance code is housed under ```thrustvtime_v4.py```. 
Both codes already have the propellants and their constants defined. However, if you wish to change the propellant analyzed, change only the constants of specific impulse, density, the a value, and the n value.
Note on the a and n value: the a and n values traditionally given need to be divided by 1000 and 10, respectively, to match the dimensions/units of the code (Ex, a value of 3.9 would need to be inputted as 0.0039, and n value of 0.40 needs to be 0.040).
The applied performance code includes the function for nozzle regression.
Both codes will return specific data in the terminal (such as burn time, max thrust, etc.) for each propellant.

# Libraries used
The ideal performance code uses the ```matplotlib``` and ```numpy``` libraries.
The applied performance code uses the same libraries as the ideal performance code, along with the ```pygasflow``` library for isentropic flow

