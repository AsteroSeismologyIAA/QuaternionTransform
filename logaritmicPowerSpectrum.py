import numpy as np
import matplotlib.pyplot as plt

c1,c2,c3 = np.loadtxt("ps.dgt",unpack=True)

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(10,10))


ax1.plot(np.log10(c1), np.log10(c2), 'b-')
ax1.set_ylabel('Log Power')

ax2.plot(np.log10(c1), np.log10(c3), 'r-')
ax2.set_xlabel('Log Frequency $\mu$Hz')
ax2.set_ylabel('Log Power')
plt.savefig("powerspectrum_18_5.png")
plt.show()
