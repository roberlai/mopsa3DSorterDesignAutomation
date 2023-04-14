import matplotlib.pyplot as plt
import numpy as np
import sys
import glob

if len(sys.argv) < 2:
    print("usage: python %s <sim_output_folder>" % sys.argv[0])
    exit(-1)

sim_output_folder = sys.argv[1]
design_name = sys.argv[1].split('/')[-3]
print(design_name)
trajectory_paths = glob.glob(sim_output_folder + '/cpp_*' +  '.txt')

startpoints = []
mean_pos = []
for dpi, trajectory_path in enumerate(trajectory_paths):

    print("trajectory_path:", trajectory_path)
    diam = 0
    with open(trajectory_path) as fp:
        fp.readline() # title
        trajectory_pos = []
        diam, has_detail_traj, num_init = \
            [float(x) for x in fp.readline().split()]
        num_init = int(num_init)

        for T in range(num_init):

            num = int(fp.readline())
            if has_detail_traj == False:
                num = 2
            data = np.zeros((int(num), 3))
            for i in range(int(num)):
                line = fp.readline()
                x, y, z = [float(v) for v in line.split()]
                data[i] = x, y, z

            if z > 0:
                continue
            trajectory_pos.append(data[-1])
            if dpi == 0:
                startpoints.append(data[0])

        trajectory_pos = np.array(trajectory_pos)
        mean_pos.append( (np.mean(trajectory_pos[:, 0]), np.mean(trajectory_pos[:, 1])) )
        print(diam, trajectory_pos[:, 0], trajectory_pos[:, 1])
        print(mean_pos[-1])
        plt.scatter(trajectory_pos[:, 0], trajectory_pos[:, 1], s=20, label="diam_" + str(int(diam)))
        plt.scatter(np.mean(trajectory_pos[:, 0]), np.mean(trajectory_pos[:, 1]), s=30, label="mean_diam_" + str(int(diam)))

startpoints = np.array(startpoints)
# plt.scatter(startpoints[:, 0], startpoints[:, 1], s=5, label="startpoint")
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
filename = design_name + "_monte_carlo.png"
plt.savefig(filename)
print('Dump to %s' % filename)
plt.show()
