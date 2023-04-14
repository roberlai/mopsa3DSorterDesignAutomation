import sys
import glob
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

DPI = 400

use_paper_font = False

if use_paper_font:
    font = {'family' : 'Times New Roman',
            'weight' : 'regular',
            'size'   : 20}
    USE_LATEX = True
    matplotlib.rc('font', **font)
    plt.rcParams.update({
        "text.usetex" : USE_LATEX,
        "font.family": "serif",
        "font.serif" : "Times"
        # "font.family" : "Helvetica",
    })

folder = sys.argv[1]

total_generation = glob.glob(folder + '/G*')

print("Total generation: %d" % len(total_generation))

top_designs = []
total_scores = []
X = []
Y = []
Gen = []
for gid in range(len(total_generation)):
    gen_path = folder + '/G' + str(gid)

    designs = glob.glob(gen_path + '/G*')

    scores = []
    score_with_name = []
    for design in designs:
        score_path = design + '/sim_output/score.txt'
        recipe_path = design + '/feature_recipe.txt'
        basename = os.path.basename(design)
        if os.path.exists(score_path):
            with open(score_path) as fp:
                score = float(fp.readline())
                scores.append(score)
                score_with_name.append( (score, basename) )

            assert(os.path.exists(recipe_path))
            with open(recipe_path) as fp:
                line = fp.readline()
                line = fp.readline()
                feature = [float(x) for x in line.strip().split()]
            Y.append([score+10])
            X.append(feature)
            Gen.append([gid])

    if len(scores) <= 0:
        break
    scores.sort(reverse=True)
    score_with_name.sort(reverse=True)
    scores = np.array(scores)
    total_scores.append(scores)
    print("Gen", gid, scores[0], np.average(scores), len(scores))
    print(score_with_name[:10])
    print(score_with_name[-10:])
    print('----------------------------------------')
    top_designs.append(score_with_name[0][1])

with open("top_designs.txt", "w") as fp:
    for name in top_designs:
        fp.write(name + '\n')

plt.xlabel("Generation")
plt.ylabel("score")
plt.boxplot(total_scores, showmeans=True, labels=[str(i) for i in range(len(total_scores))])
plt.tight_layout()
plt.savefig("ga_score.png", dpi=DPI)
