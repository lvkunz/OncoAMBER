from Process import *
from scipy.stats import beta

r = beta.rvs(1.5, 1.5, size=1000)

plt.hist(r, bins=20)
plt.show()
