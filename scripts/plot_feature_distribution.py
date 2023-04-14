import sys
import glob
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

DPI = 400
SHOW_MEANS = True

use_paper_font = False

if use_paper_font:
    font = {'family' : 'Times New Roman',
           'weight' : 'regular',
           'size'   : 20}
    matplotlib.rc('font', **font)
    plt.rcParams.update({
       "text.usetex" : True,
       "font.family": "serif",
       "font.serif" : "Times",
       "font.family" : "Helvetica",
    })

folder = sys.argv[1]
gen_id = int(sys.argv[2])
read_score = int(sys.argv[3])

total_desgins = glob.glob(folder + '/G' + str(gen_id) + '/*')

print("Size of generation %d: %d" % (gen_id, len(total_desgins)))

features = [ [] for i in range(7) ]
# feature_name = []
for design in total_desgins:
    recipe = design + '/feature_recipe.txt'
    with open(recipe) as fp:
        _ = fp.readline().strip().split(',')
        feature = [float(x) for x in fp.readline().split() ]
        feature.pop(6)
        feature.pop(7)

    for i in range(7):
        features[i].append( feature[i] )

feature_name = ["" for _ in range(7)]
feature_name[0] = '$radius$'
feature_name[1] = '$theta$'
feature_name[2] = '$phi$'
feature_name[3] = '$y_{init}$'
feature_name[4] = '$y_{init\_shift}$'
feature_name[5] = '$y_{gap}$'
feature_name[6] = '$z_{gap}$'

title = "$G" + str(gen_id) + "$"

features = np.array(features)
plt.title(title)
plt.boxplot(list(features), showmeans=SHOW_MEANS, labels=feature_name)
plt.gcf().autofmt_xdate()
plt.savefig("feature_dis_G%d.png" % (gen_id), dpi=DPI)

if read_score:
    scores = []
    for design in total_desgins:

        score_path = design + '/sim_output/score.txt'
        score = 0
        with open(score_path) as fp:
            score = float(fp.readline())
        scores.append(score)

    scores = np.array(scores)

    top = np.argsort(scores)[-100:][::-1]
    print(top)
    print(scores[top])
    print(list(features[:, top]))

    plt.clf()
    plt.title(title)
    plt.boxplot(list(features[:, top]), showmeans=SHOW_MEANS, labels=feature_name)
    plt.gcf().autofmt_xdate()
    plt.savefig("feature_dis_G%d_top100.png" % (gen_id), dpi=DPI)
