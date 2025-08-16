import numpy as np
import matplotlib.pyplot as plt
print("Thrust vs Time Simulation")
# Constants
g= 9.80665

# Propellant constants
# Propellant 1: AP 70% Al 10% HTPB 17% PbSt 3%
# Propellant 2: AP 65% Al 15% HTPB 15% PbSt 5%
# Propellant 3: AP 70% Al 10% HTPB 17% FeO 3%
# Propellant 4: AP 68% Al 15% HTPB 12% FeO 5%

specific_impulse_1 = 250
density_1 = 1770
a_const_1 = 0.0036
n_const_1 = 0.032

specific_impulse_2 = 245
density_2 = 1780
a_const_2 = 0.0040
n_const_2 = 0.030

specific_impulse_3 = 250
density_3 = 1740
a_const_3 = 0.0032
n_const_3 = 0.038

specific_impulse_4 = 255
density_4 = 1800
a_const_4 = 0.0039
n_const_4 = 0.038

exit_velocity_1 = specific_impulse_1 * g
exit_velocity_2 = specific_impulse_2 * g
exit_velocity_3 = specific_impulse_3 * g
exit_velocity_4 = specific_impulse_4 * g

x_time_1 = []
x_time_2 = []
x_time_3 = []
x_time_4 = []

y_thrust_1 = []
y_thrust_2 = []
y_thrust_3 = []
y_thrust_4 = []


# time + rad steps, do time increments on graph / 10 at least, ex 0.01 increment = 0.001 time step
time_step = 0.1
radius_step = 0.1

# Constants for burn area / Grain geometry
grain_length = 19
r_initial = 0.05
r_final = 0.809
r_leeway = 0.3
sqrt_constant = (r_final-r_leeway-r_initial) / np.sqrt(grain_length)

#Constants for burn rate
throat_area = 1
# 0.1842434 <- Northrop Grumman's throat area w/ area ratio of 10 assumed
burn_area_initial = 127.05294011645312
c_star = 1500

chamber_pressure_initial_1 = np.power(((a_const_1 * density_1 * burn_area_initial * c_star) / throat_area), (1/(1-n_const_1)))
chamber_pressure_initial_2 = np.power(((a_const_2 * density_2 * burn_area_initial * c_star) / throat_area), (1/(1-n_const_2)))
chamber_pressure_initial_3 = np.power(((a_const_3 * density_3 * burn_area_initial * c_star) / throat_area), (1/(1-n_const_3)))
chamber_pressure_initial_4 = np.power(((a_const_4 * density_4 * burn_area_initial * c_star) / throat_area), (1/(1-n_const_4)))

burn_rate_initial_1 = a_const_1 * np.power(chamber_pressure_initial_1, n_const_1)
burn_rate_initial_2 = a_const_2 * np.power(chamber_pressure_initial_2, n_const_2)
burn_rate_initial_3 = a_const_3 * np.power(chamber_pressure_initial_3, n_const_3)
burn_rate_initial_4 = a_const_4 * np.power(chamber_pressure_initial_4, n_const_4)

r_array = []

burn_area_array = []
burn_area_array_1 = []
burn_area_array_2 = []
burn_area_array_3 = []
burn_area_array_4 = []

impulse_array = []
impulse_array_1 = []
impulse_array_2 = []
impulse_array_3 = []
impulse_array_4 = []


for j in np.arange(0, grain_length, radius_step):
    r_array.append((sqrt_constant)*np.sqrt(j)+ r_initial)
    #r_array.append((0.03182*j)+ r_initial) # linear grain geometry






# PROPELLANT 1 CALCULATIONS
r_array_1 = r_array.copy()
t=0
while(r_array_1[0] < r_final):
    for i in range(len(r_array_1)):
        if(t!=0):
            r_array_1[i] = r_array_1[i] + burn_rate_1 * time_step
        if (t == 0):
            burn_rate_1 = burn_rate_initial_1
        if r_array_1[i] > r_final:
            r_array_1[i] = r_final
        if r_array_1[i] < 0:
            r_array_1[i] = 0
        # Cylindrical grain surface area
        #burn_area = 2 * np.pi * r_array_1[i] * radius_step

        # Star grain approximation
        epsilon = 2.0  # extra area enhancement factor (e.g., 100% more at t=0)
        inst_burn_time = (r_final - r_array_1[i]) / burn_rate_1  # instantaneous burn time
        if((r_final - r_array_1[i]) == 0):
            inst_burn_time = 1
        beta = 1 / inst_burn_time  # decay rate of star shape
        eta_t = 1 + epsilon * np.exp(-beta * t)
        burn_area_1 = 2 * np.pi * r_array_1[i] * radius_step * eta_t
        burn_area_array_1.append(burn_area_1)

        if r_array_1[i] >= r_final:
            burn_area_1 = 0

    #prep for next time step
    chamber_pressure_1 = np.power((a_const_1 * density_1 * sum(burn_area_array_1)*c_star), (1/(1-n_const_1)))
    burn_rate_1 = a_const_1 * np.power(chamber_pressure_1, n_const_1)
    mass_flow_rate_1 = sum(burn_area_array_1) * density_1 * burn_rate_1
    thrust_1 = mass_flow_rate_1 * exit_velocity_1
    if(t==0):
        print("Propellant 1 Initial burn area: " + str(sum(burn_area_array_1)))
    #supposed_throat_area = (sum(burn_area_array_1) * density_1 * burn_rate_1 * c_star) / chamber_pressure_1
    #print(supposed_throat_area)
    burn_area_array_1.clear()
    y_thrust_1.append(thrust_1)
    x_time_1.append(t)
    impulse_array_1.append(thrust_1 * time_step)
    t += time_step






# PROPELLANT 2 CALCULATIONS
r_array_2 = r_array.copy()
t=0
while(r_array_2[0] < r_final):
    for i in range(len(r_array_2)):
        if(t!=0):
            r_array_2[i] = r_array_2[i] + burn_rate_2 * time_step
        if (t == 0):
            burn_rate_2 = burn_rate_initial_2
        if r_array_2[i] > r_final:
            r_array_2[i] = r_final
        if r_array_2[i] < 0:
            r_array_2[i] = 0
        # Cylindrical grain surface area
        #burn_area = 2 * np.pi * r_array_2[i] * radius_step

        # Star grain approximation
        epsilon = 2.0  # extra area enhancement factor (e.g., 100% more at t=0)
        inst_burn_time = (r_final - r_array_2[i]) / burn_rate_2  # instantaneous burn time
        if((r_final - r_array_2[i]) == 0):
            inst_burn_time = 1
        beta = 1 / inst_burn_time  # decay rate of star shape
        eta_t = 1 + epsilon * np.exp(-beta * t)
        burn_area_2 = 2 * np.pi * r_array_2[i] * radius_step * eta_t
        burn_area_array_2.append(burn_area_2)

        if r_array_2[i] >= r_final:
            burn_area_2 = 0

    #prep for next time step
    chamber_pressure_2 = np.power((a_const_2 * density_2 * sum(burn_area_array_2)*c_star), (1/(1-n_const_2)))
    burn_rate_2 = a_const_2 * np.power(chamber_pressure_2, n_const_2)
    mass_flow_rate_2 = sum(burn_area_array_2) * density_2 * burn_rate_2
    thrust_2 = mass_flow_rate_2 * exit_velocity_2
    #supposed_throat_area = (sum(burn_area_array_2) * density_2 * burn_rate_2 * c_star) / chamber_pressure_2
    #print(supposed_throat_area)
    if(t==0):
        print("Propellant 2 Initial burn area: " + str(sum(burn_area_array_2)))
    burn_area_array_2.clear()
    y_thrust_2.append(thrust_2)
    x_time_2.append(t)
    impulse_array_2.append(thrust_2 * time_step)
    t += time_step






# PROPELLANT 3 CALCULATIONS
r_array_3 = r_array.copy()
t=0
while(r_array_3[0] < r_final):
    for i in range(len(r_array_3)):
        if(t!=0):
            r_array_3[i] = r_array_3[i] + burn_rate_3 * time_step
        if (t == 0):
            burn_rate_3 = burn_rate_initial_3
        if r_array_3[i] > r_final:
            r_array_3[i] = r_final
        if r_array_3[i] < 0:
            r_array_3[i] = 0
        # Cylindrical grain surface area
        #burn_area = 2 * np.pi * r_array_3[i] * radius_step

        # Star grain approximation
        epsilon = 2.0  # extra area enhancement factor (e.g., 100% more at t=0)
        inst_burn_time = (r_final - r_array_3[i]) / burn_rate_3  # instantaneous burn time
        if((r_final - r_array_3[i]) == 0):
            inst_burn_time = 1
        beta = 1 / inst_burn_time  # decay rate of star shape
        eta_t = 1 + epsilon * np.exp(-beta * t)
        burn_area_3 = 2 * np.pi * r_array_3[i] * radius_step * eta_t
        burn_area_array_3.append(burn_area_3)

        if r_array_3[i] >= r_final:
            burn_area_3 = 0

    #prep for next time step
    chamber_pressure_3 = np.power((a_const_3 * density_3 * sum(burn_area_array_3)*c_star), (1/(1-n_const_3)))
    burn_rate_3 = a_const_3 * np.power(chamber_pressure_3, n_const_3)
    mass_flow_rate_3 = sum(burn_area_array_3) * density_3 * burn_rate_3
    thrust_3 = mass_flow_rate_3 * exit_velocity_3
    #supposed_throat_area = (sum(burn_area_array_3) * density_3 * burn_rate_3 * c_star) / chamber_pressure_3
    #print(supposed_throat_area)
    if(t==0):
        print("Propellant 3 Initial burn area: " + str(sum(burn_area_array_3)))
    burn_area_array_3.clear()
    y_thrust_3.append(thrust_3)
    x_time_3.append(t)
    impulse_array_3.append(thrust_3 * time_step)
    t += time_step






# PROPELLANT 4 CALCULATIONS
r_array_4 = r_array.copy()
t=0
while(r_array_4[0] < r_final):
    for i in range(len(r_array_4)):
        if(t!=0):
            r_array_4[i] = r_array_4[i] + burn_rate_4 * time_step
        if (t == 0):
            burn_rate_4 = burn_rate_initial_4
        if r_array_4[i] > r_final:
            r_array_4[i] = r_final
        if r_array_4[i] < 0:
            r_array_4[i] = 0
        # Cylindrical grain surface area
        #burn_area = 2 * np.pi * r_array_4[i] * radius_step

        # Star grain approximation
        epsilon = 2.0  # extra area enhancement factor
        inst_burn_time = (r_final - r_array_4[i]) / burn_rate_4  # instantaneous burn time
        if((r_final - r_array_4[i]) == 0):
            inst_burn_time = 1
        beta = 1 / inst_burn_time  # decay rate of star shape
        eta_t = 1 + epsilon * np.exp(-beta * t)
        burn_area_4 = 2 * np.pi * r_array_4[i] * radius_step * eta_t
        burn_area_array_4.append(burn_area_4)

        if r_array_4[i] >= r_final:
            burn_area_4 = 0

    #prep for next time step
    chamber_pressure_4 = np.power((a_const_4 * density_4 * sum(burn_area_array_4)*c_star), (1/(1-n_const_4)))
    burn_rate_4 = a_const_4 * np.power(chamber_pressure_4, n_const_4)
    mass_flow_rate_4 = sum(burn_area_array_4) * density_4 * burn_rate_4
    thrust_4 = mass_flow_rate_4 * exit_velocity_4
    #supposed_throat_area = (sum(burn_area_array_4) * density_4 * burn_rate_4 * c_star) / chamber_pressure_4
    #print(supposed_throat_area)
    if(t==0):
        print("Propellant 4 Initial burn area: " + str(sum(burn_area_array_4)))
    burn_area_array_4.clear()
    y_thrust_4.append(thrust_4)
    x_time_4.append(t)
    impulse_array_4.append(thrust_4 * time_step)
    t += time_step


print("Propellant 1 Initial Chamber pressure: " + str(chamber_pressure_initial_1))
print("Propellant 2 Initial Chamber pressure: " + str(chamber_pressure_initial_2))
print("Propellant 3 Initial Chamber pressure: " + str(chamber_pressure_initial_3))
print("Propellant 4 Initial Chamber pressure: " + str(chamber_pressure_initial_4))

print("Propellant 1 Initial Burn Rate: " + str(burn_rate_initial_1))
print("Propellant 2 Initial Burn Rate: " + str(burn_rate_initial_2))
print("Propellant 3 Initial Burn Rate: " + str(burn_rate_initial_3))
print("Propellant 4 Initial Burn Rate: " + str(burn_rate_initial_4))

print("Propellant 1 Burn Time: " + str(len(x_time_1) * time_step))
print("Propellant 2 Burn Time: " + str(len(x_time_2) * time_step))
print("Propellant 3 Burn Time: " + str(len(x_time_3) * time_step))
print("Propellant 4 Burn Time: " + str(len(x_time_4) * time_step))

print("Propellant 1 Max Thrust: " + str(max(y_thrust_1)))
print("Propellant 2 Max Thrust: " + str(max(y_thrust_2)))
print("Propellant 3 Max Thrust: " + str(max(y_thrust_3)))
print("Propellant 4 Max Thrust: " + str(max(y_thrust_4)))

print("Propellant 1 Total Impulse: " + str(sum(impulse_array_1)))
print("Propellant 2 Total Impulse: " + str(sum(impulse_array_2)))
print("Propellant 3 Total Impulse: " + str(sum(impulse_array_3)))
print("Propellant 4 Total Impulse: " + str(sum(impulse_array_4)))

print("Propellant 1 Average Thrust: " + str(sum(y_thrust_1) / len(y_thrust_1)))
print("Propellant 2 Average Thrust: " + str(sum(y_thrust_2) / len(y_thrust_2)))
print("Propellant 3 Average Thrust: " + str(sum(y_thrust_3) / len(y_thrust_3)))
print("Propellant 4 Average Thrust: " + str(sum(y_thrust_4) / len(y_thrust_4)))

plt.figure()
plt.plot(x_time_1, y_thrust_1, label='Propellant 1')
plt.plot(x_time_2, y_thrust_2, label='Propellant 2')
plt.plot(x_time_3, y_thrust_3, label='Propellant 3')
plt.plot(x_time_4, y_thrust_4, label='Propellant 4')
plt.title("Thrust vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Thrust (N)")
plt.legend()
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.grid()
plt.show()
