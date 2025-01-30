import numpy as np
import matplotlib.pyplot as plt
import sys

PowerDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/GetFluence/Data/Power/"
SavePath = "/Users/chiche/Desktop/PowerMaps/"
ZenithAll, EnergyAll, DepthAll, EfieldAir_x, EfieldAir_y, EfieldAir_z, EfieldAir_tot, fAir_x, fAir_y, fAir_z, fAir_tot = np.loadtxt(PowerDataPath + "AirPowerMap.txt", unpack =True)

ZenithAll, Energy, Depth, EfieldIce_x, EfieldIce_y, EfieldIce_z, EfieldIce_tot, fIce_x, fIce_y, fIce_z, fIce_tot = np.loadtxt(PowerDataPath + "IcePowerMap.txt", unpack = True)
EnergyAll = np.array([0.0316, 0.1, 0.316])
glevel =3216
DepthAll = glevel - np.array([0,40, 60, 80, 100])


for i in range(len(EnergyAll)):
    for j in range(len(DepthAll)):
        
        sel = ((Energy== EnergyAll[i]) & (Depth== DepthAll[j]))

        #In-air emission
        plt.plot(ZenithAll[sel], EfieldAir_x[sel], label ="Ex")
        plt.plot(ZenithAll[sel], EfieldAir_y[sel], label ="Ey")
        plt.plot(ZenithAll[sel], EfieldAir_z[sel], label ="Ez")
        #plt.plot(ZenithAll[sel], EfieldAir_tot[sel], label ="Etot")
        plt.xlabel("zenith [Deg.]")
        plt.ylabel("Integrated Efield $\mu Vs/m$")
        plt.title("In-air, E=%.2f, |z|=%.d" %(EnergyAll[i], DepthAll[j]))
        plt.legend()
        plt.yscale("log")
        plt.savefig(SavePath + "AirEfieldMap_E%.2f_z%.d.pdf" %(EnergyAll[i], DepthAll[j]), bbox_inches ="tight")
        plt.show()
        plt.close()

        
        plt.plot(ZenithAll[sel], fAir_x[sel], label ="Power-x")
        plt.plot(ZenithAll[sel], fAir_y[sel], label ="Power-y")
        plt.plot(ZenithAll[sel], fAir_z[sel], label ="Power-z")
        #plt.plot(ZenithAll[sel], fAir_tot[sel], label ="Etot")
        plt.xlabel("zenith [Deg.]")
        plt.ylabel("Integrated fluence $(\mu V)^{2}\, s$")
        plt.title("In-air, E=%.2f, |z|=%.d" %(EnergyAll[i], DepthAll[j]))
        plt.legend()
        plt.yscale("log")
        plt.savefig(SavePath + "AirPowerMap_E%.2f_z%.d.pdf" %(EnergyAll[i], DepthAll[j]))
        plt.show()
        plt.close()
        

        #In-ice emission
        plt.plot(ZenithAll[sel], EfieldIce_x[sel], label ="Ex")
        plt.plot(ZenithAll[sel], EfieldIce_y[sel], label ="Ey")
        plt.plot(ZenithAll[sel], EfieldIce_z[sel], label ="Ez")
        #plt.plot(ZenithAll[sel], EfieldIce_tot[sel], label ="Etot")
        plt.xlabel("zenith [Deg.]")
        plt.ylabel("Integrated Efield $\mu Vs/m$")
        plt.title("In-ice, E=%.2f, |z|=%.d" %(EnergyAll[i], DepthAll[j]))
        plt.legend()
        plt.yscale("log")
        plt.savefig(SavePath + "IceEfieldMap_E%.2f_z%.d.pdf" %(EnergyAll[i], DepthAll[j]), bbox_inches ="tight")
        plt.show()
        plt.close()

        plt.plot(ZenithAll[sel], fIce_x[sel], label ="Power-x")
        plt.plot(ZenithAll[sel], fIce_y[sel], label ="Power-y")
        plt.plot(ZenithAll[sel], fIce_z[sel], label ="Power-z")
        #plt.plot(ZenithAll[sel], fIce_tot[sel], label ="Etot")
        plt.xlabel("zenith [Deg.]")
        plt.ylabel("Integrated fluence $(\mu V)^{2}\, s$")
        plt.title("In-ice, E=%.2f, |z|=%.d" %(EnergyAll[i], DepthAll[j]))
        plt.legend()
        plt.yscale("log")
        plt.savefig(SavePath + "IcePowerMap_E%.2f_z%.d.pdf" %(EnergyAll[i], DepthAll[j]))
        plt.show()
        plt.close()

        #sys.exit()

























#--------------------------------------------#
"""
PowerDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/GetFluence/Data/Power/"
zenith, AirPow, IcePow = np.loadtxt(PowerDataPath + "power.txt", unpack =True)


plt.plot(zenith, AirPow, label = "In-air power 100m depth")
plt.plot(zenith, IcePow, label = "In-ice power 100m depth")
plt.legend()
plt.xlabel("zenith [Deg]")
plt.show()

plt.scatter(zenith, IcePow/AirPow, label = "In-air power 100m depth")
plt.legend()
plt.xlabel("zenith [Deg]")
plt.ylabel("Ice/Air power")
plt.yscale("log")
plt.show()


zenith, AirPow, IcePow = np.loadtxt(PowerDataPath + "fluencevs_theta.txt", unpack =True)


plt.plot(zenith, AirPow, label = "In-air power 100m depth")
plt.plot(zenith, IcePow, label = "In-ice power 100m depth")
plt.legend()
plt.xlabel("zenith [Deg]")
plt.show()

plt.scatter(zenith, IcePow/AirPow, label = "In-air power 100m depth")
plt.legend()
plt.xlabel("zenith [Deg]")
plt.ylabel("Ice/Air power")
plt.yscale("log")
plt.show()



zenith, AirPow, IcePow = np.loadtxt(PowerDataPath + "Efield_vs_theta.txt", unpack =True)


plt.plot(zenith, AirPow, label = "In-air power 100m depth")
plt.plot(zenith, IcePow, label = "In-ice power 100m depth")
plt.legend()
plt.xlabel("zenith [Deg]")
plt.show()

plt.scatter(zenith, IcePow/AirPow, label = "In-air power 100m depth")
plt.legend()
plt.xlabel("zenith [Deg]")
plt.ylabel("Ice/Air power")
plt.yscale("log")
plt.show()

### plot for different polarizations
"""