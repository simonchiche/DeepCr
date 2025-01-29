import numpy as np
import matplotlib.pyplot as plt

#zenith, AirPow, IcePow = np.loadtxt("./fluencevs_theta.txt", unpack =True)
#zenith, AirPow, IcePow = np.loadtxt("./Efield_vs_theta.txt", unpack =True)
zenith, AirPow, IcePow = np.loadtxt("./power.txt", unpack =True)


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


zenith, AirPow, IcePow = np.loadtxt("./fluencevs_theta.txt", unpack =True)


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



zenith, AirPow, IcePow = np.loadtxt("./Efield_vs_theta.txt", unpack =True)


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