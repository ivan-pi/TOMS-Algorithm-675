import numpy as np
import matplotlib.pyplot as plt

temperature = np.asarray([15.5,28.5,42,54,62,77,66,66,66,13])
relhumid = np.asarray([16.5,24,24,34,35.5,38,69.5,43.5,69.5,64])
time = np.asarray([0,7.5,3,3.5,3,22.5,23,5.5,59.5,5.5,3.5])

time_factor = 157.5/5.
time /= time_factor
time = np.cumsum(time)

rf_factor = 46.5/30.
rf0 = 100 - 79.5/rf_factor

relhumid = rf0 + relhumid/rf_factor

temp_factor = 62./40.
temp0 = 90 - 79.5/temp_factor

temperature = temp0 + temperature/temp_factor

def rf(t):
    conds = [t >= tt for tt in time[:-1]]
    return np.piecewise(t,conds,relhumid)

def temp(t):
    conds = [t >= tt for tt in time[:-1]]
    return np.piecewise(t,conds,temperature)


def main():

    t = np.linspace(0,5,400)
    r = rf(t)
    tt = temp(t)

    print("Times = ", time)
    print("Relative humidity = ", np.round(relhumid))
    print("temperature = ", np.round(temperature))

    plt.plot(t,r)
    plt.plot(time[0:-1],relhumid,'o')

    plt.plot(t,tt)
    plt.plot(time[0:-1],temperature,'o')
    plt.show()

if __name__ == '__main__':
    main()