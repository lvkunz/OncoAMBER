import matplotlib.pyplot as plt

from Vessel import *

vasculature = VasculatureNetwork()
vasculature.build_vasculature(Sphere(center = np.array([0,0,0]), radius = 3.0).generate_random_points(1000))
vasculature.save('Vasculature/vasculature_1000.txt')

# vasculature = VasculatureNetwork()
# vasculature.read('Vasculature/vasculature.txt')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
vasculature.plot(fig, ax)
plt.savefig('Vasculature/vasculature.png')
plt.show()

