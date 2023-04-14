import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import sys
import glob
import numpy as np

from vispy import app, visuals, scene
from vispy.visuals.transforms import STTransform
from vispy.visuals.markers import MarkersVisual

if len(sys.argv) < 2:
    print("usage: python %s <sim_output_folder>"
        " [num_traj=100]" % sys.argv[0]
    )
    exit(-1)

sim_output_folder = sys.argv[1]
print("sim_outut_folder:", sim_output_folder)

obstacle_ndoes_paths = glob.glob(sim_output_folder + '/*obstacle_nodes.log')

if len(obstacle_ndoes_paths):
    obstacle_ndoes_paths = obstacle_ndoes_paths[0]
    print("obstacle_ndoes_paths:", obstacle_ndoes_paths)
else:
    print("Cannot find obstacle ndoes log in %s" % sim_output_folder)
    exit(-1)

num_traj = 100
if len(sys.argv) >= 3:
    num_traj = int(sys.argv[2])

print("num_traj:", num_traj)
        
obstacle_pos = []
with open(obstacle_ndoes_paths) as fp:
    num = int(fp.readline())
    obstacle_pos = np.zeros((num,3))
    for i in range(num):
        x, y, z = [float(d) for d in fp.readline().split()]
        obstacle_pos[i] = x, y, z
        
### Start Plotting
Scatter3D = scene.visuals.create_visual_node(visuals.MarkersVisual)

canvas = scene.SceneCanvas(keys="interactive", show=True)

view = canvas.central_widget.add_view()
view.camera = "turntable"
view.camera.fov = 45
view.camera.distance = 500

p1 = Scatter3D(parent=view.scene)
p1.set_gl_state("translucent", blend=True, depth_test=True)

p1.set_data(
    obstacle_pos, symbol="o", size=1, edge_width=0.5, edge_color="blue"
)

colors = ['orange', 'green', 'red', 'black', 'purple', 'yellow']
trajectory_paths = glob.glob(sim_output_folder + '/cpp_*' +  '.txt')
target_T = 0
for dpi, trajectory_path in enumerate(trajectory_paths):

    print("trajectory_path:", trajectory_path)
    trajectory_pos = []
    diam = 0
    with open(trajectory_path) as fp:
        fp.readline() # title
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
            trajectory_pos.append(data)

    if len(trajectory_pos[target_T]) < num_traj:
        step = 1
    else: 
        step = int(len(trajectory_pos[target_T]) / num_traj)

    for i in range(0, len(trajectory_pos[target_T]), step):
        sphere1 = scene.visuals.Sphere(
            radius=diam/2.0, method='ico', parent=view.scene,
            edge_color=colors[dpi]
        )
        sphere1.transform = STTransform(translate=trajectory_pos[target_T][i])

if __name__ == "__main__":
    if sys.flags.interactive != 1:
        app.run()
