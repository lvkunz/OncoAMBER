import matplotlib.pyplot as plt

from Vessel import *

vasculature = VasculatureNetwork()
vasculature.generate_vasculature(5000)
vasculature.save('Vasculature/vasculature_5000.txt')

# vasculature = VasculatureNetwork()
# vasculature.read('Vasculature/vasculature.txt')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
vasculature.plot(fig, ax)
# plt.savefig('Vasculature/vasculature.png')
plt.show()

