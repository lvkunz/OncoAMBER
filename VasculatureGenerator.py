import matplotlib.pyplot as plt

from Vessel import *

vasculature = Vasculature()
vasculature.generate_vasculature(2000)
vasculature.save('Vasculature/vasculature.txt')

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# vasculature.plot(fig, ax)
# plt.savefig('vasculature.png')
# plt.show()

