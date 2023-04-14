import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np

if len(sys.argv) < 2:
    print("usage: python %s <file>" % sys.argv[0])
    exit(-1)

mpl.rcParams['figure.dpi'] = 300

xs = []
ys = []
zs = []
with open(sys.argv[1]) as fp:
    num = int(fp.readline())
    for i in range(num):
        x, y, z = [float(d) for d in fp.readline().split()]
        xs.append(x)
        ys.append(y)
        zs.append(z)
        
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.view_init(23, 23)

ax.scatter(xs, ys, zs, s=0.5)

X = np.array(xs)
Y = np.array(ys)
Z = np.array(zs)

max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())

for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')

class sphere:

    def __init__(self, radius, x, y, z):

        d = 40
        theta = np.linspace(0, 2 * np.pi, d)
        phi   = np.linspace(0, np.pi, d/2)

        self.xs = x + radius * np.outer(np.sin(phi), np.cos(theta))
        self.ys = y + radius * np.outer(np.sin(phi), np.sin(theta))
        self.zs = z + radius * np.outer(np.cos(phi), np.ones(theta.size))

    def draw(self, ax):
        # ax.scatter(
        #   self.xs.flatten(), 
        #   self.ys.flatten(), 
        #   self.zs.flatten(), s= 0.5)
        return ax.plot_surface(self.xs, self.ys, self.zs, color='red')
        # ax.plot_wireframe(self.xs, self.ys, self.zs, color ='red')
        
diam = 15
#sphere(diam/2.0, 45, 0, 75).draw(ax)
plt.show()
