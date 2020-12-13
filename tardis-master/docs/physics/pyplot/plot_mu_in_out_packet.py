from pylab import *
from astropy import units as u, constants as const

x, y = x, y = mgrid[1:1000, 1:1000]
mu_in = x / 500.0 - 1
mu_out = y / 500.0 - 1
v = 1.1e4 * u.km / u.s
doppler_fac = (1 - mu_in * v / const.c) / (1 - mu_out * v / const.c)
xlabel("$\mu_{\\rm in}$")
ylabel("$\mu_{\\rm out}$")
imshow(np.rot90(doppler_fac), extent=[-1, 1, -1, 1], cmap="bwr")
colorbar()
show()
