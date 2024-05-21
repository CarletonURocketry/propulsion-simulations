import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize
from CoolProp.CoolProp import PropsSI
#===========================| INITIAL CONDITIONS |===================================================================

T_0 = 3+273 #K
M_intial = 6.3 #kg. initial amount of nitrous loaded. set to 0 to have this solved for you as a function of a given quality and tempreature
x_initial = 0.1 #initial quality of nitrous. ignore if using m_initial
gas = 'CO2' #either CO2 or N2O
V_tank = 0.012#1.25 * math.pi / 4 * (0.15) ** 2  # m^3
Cd_injector = 0.08# Discharge Coefficient of individual injector elements. see pg 279 of rocket propulsion elements
A_injector = 0.000106029#(1.5*10**-3)**2/4*math.pi*60
P_chamber = 101325#450/14.5*101325 #should be the pressure directly downstream of the injector
timespan = 10#s
timestep = 0.001#s

#==========================================| FUNCTIONS |===========================================================
def Temperature(T,U_total,M_total):
    #this function takes a guess for the temperature (T) and the total energy and mass in the tank and returns a value near 0 when the temperature is correct
    if T >=309:
        T=309
    elif T<= 183:
        T=183
    u_liq = PropsSI('U','T',T,'Q',0,gas)
    u_vap = PropsSI('U','T',T,'Q',1,gas)
    rho_liq = PropsSI('D','T',T,'Q',0,gas)
    rho_gas = PropsSI('D','T',T,'Q',1,gas)
    x = ((U_total/M_total) - u_liq)/(u_vap-u_liq) # find the quality based off the internal energy of the tank
    Zero = V_tank - M_total*((1-x)/rho_liq+x/rho_gas) # use this equality and density to find the temperature that satisifies this eqn
    return Zero
#def Z_solver(T,m):

def N2O_Properties (T,M_total):
    #this function takes the temperature of the nitrous and the total mass and returns various properties of the nitrous in the tank
    rho_liq = PropsSI('D', 'T', T, 'Q', 0, gas)
    rho_gas = PropsSI('D', 'T', T, 'Q', 1, gas)
    h_liq = PropsSI('H','T',T,'Q',0,gas)
    P = PropsSI('P','T',T,'Q',1,gas)
    x = ((V_tank/M_total)-(1/rho_liq))/((1/rho_gas)-(1/rho_liq))
    m_vap = x*M_total
    m_liq = (1-x)*M_total
    return rho_liq,rho_gas,h_liq,P,x,m_liq,m_vap

def Differential_Equation(time,Values):
    #this equation yeilds values for the rate of change of mass and energy in the tank by taking the current mass and energy in the tank, solving for the temperature and then using the liquid properties from the temperature to find the mass flow and energy flow
    #Values is a tuple containing the total energy and total mass, and t is the current time, this fucntion gets integrated with respect to time to get a function for the tank mass and energy over time
    U_total, M_total = Values
    T= scipy.optimize.root(Temperature,x0=280,args=(U_total,M_total),method='lm', tol=1e-7).x
    rho_liq,rho_gas,h_liq,P,x,m_liq,m_vap = N2O_Properties(T,M_total)
    if m_liq <=0.1:
        m_dot = 0
        #when the liquid mass is near zero, the mass flow is set to zero and the burn is considered "done". a gas-phase sim should be added in place of this later on
    else:
        m_dot = Cd_injector*A_injector*math.sqrt(2*rho_liq*(P-P_chamber))
    Mass_dot_total = -m_dot
    U_dot_total = -m_dot*h_liq
    return  (U_dot_total,Mass_dot_total)
#create an array of times from 0 to the chosen duration with the chosen timestep
time_list = np.linspace(0,timespan,int(timespan/timestep))
#create empty arrays to assign numbers from the solution
m_liq_list = np.zeros(len(time_list))
m_vap_list = np.zeros(len(time_list))
P_list = np.zeros(len(time_list))
rho_liq_list=np.zeros(len(time_list))
rho_gas_list=np.zeros(len(time_list))
x_list = np.zeros(len(time_list))
h_liq_list = np.zeros(len(time_list))
m_dot_list = np.zeros(len(time_list))
Temperature_list = np.zeros(len(time_list))
#=================================|FINDING TANK MASS|======================================

#================================|SOLVING THE EQUATION|=========================================
rho_liq_0 , rho_gas_0 , h_liq_0, P_0,x_0,m_liq_0,m_vap_0 = N2O_Properties(T_0,M_intial)
if x_0 > 1:
    print('too much mass in tank')
print(rho_liq_0)
U_initial = M_intial*PropsSI('U','T',T_0,'Q',x_0,gas)
values_initial = (U_initial,M_intial)
the_answer = scipy.integrate.odeint(Differential_Equation,y0=values_initial,t=time_list,tfirst=True)
for i in range(0,len(the_answer)):
    U_total, M_total =the_answer.T[0][i],the_answer.T[1][i]
    T= scipy.optimize.root(Temperature,x0=280,args=(U_total,M_total),method='lm', tol=1e-7).x
    rho_liq_list[i], rho_gas_list[i],h_liq_list[i],P_list[i],x_list[i],m_liq_list[i],m_vap_list[i] = N2O_Properties(T,M_total)
    m_dot_list[i] = A_injector*Cd_injector*math.sqrt(rho_liq_list[i]*2*(P_list[i]-P_chamber))
    Temperature_list[i] = T-273
plt.plot(time_list, Temperature_list)
print(m_dot_list)
# Adding labels and title
plt.xlabel('Time (seconds)')
plt.ylabel('Temperature (C)')
plt.title('Tank Temperature vs Time')

# Show the plot
plt.show()
plt.plot(time_list, P_list)

# Adding labels and title
plt.xlabel('Time (seconds)')
plt.ylabel('Pressure (Pa)')
plt.title('Tank Pressure vs Time')

# Show the plot
plt.show()

plt.plot(time_list, P_list)

# Adding labels and title
plt.xlabel('Time (seconds)')
plt.ylabel('Pressure (Pa)')
plt.title('Tank Pressure vs Time')

# Show the plot
plt.show()
plt.plot(time_list, m_liq_list, label='Liquid Mass')

# Plot m_vap_list against time_list
plt.plot(time_list, m_vap_list, label='Vapor Mass')

# Set labels and title
plt.xlabel('Time')
plt.ylabel('Mass')
plt.title('Liquid and Vapor Mass vs. Time')

# Add legend
plt.legend()

# Show plot
plt.show()
plt.plot(time_list, m_dot_list)

# Adding labels and title
plt.xlabel('Time (seconds)')
plt.ylabel('Mass Flow (kg/s)')
plt.title('Tank Mass Flow vs Time')

# Show the plot
plt.show()
