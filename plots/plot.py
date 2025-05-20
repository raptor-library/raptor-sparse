import numpy as np
import glob

computer = "dane"
leaders = [4, 10, 20]

import pyfancyplot as plt
import matplotlib
matplotlib.use("qtagg")

class SpMV:
    def __init__(self):
        self.time = np.inf

    def add_time(self, time):
        if time < self.time:
            self.time = time


class Matrix:
    name = ""
    
    def __init__(self, name):
        self.name = name
        self.std = list()
        self.early_std = list()
        self.csc = list()
        self.early_csc = list()
        self.n = 0

    def add_test(self, n):
        if (self.n != n):
            self.n = n
            self.std.append(SpMV())
            self.early_std.append(SpMV())
            self.csc.append(SpMV())
            self.early_csc.append(SpMV())

    def add_time(self, name, time):
        if "Standard" in name:
            self.std[-1].add_time(time)
        elif "Early Comp" in name:
            self.early_std[-1].add_time(time)
        elif "Early Recv" in name:
            self.early_csc[-1].add_time(time)
        elif "CSC" in name:
            self.csc[-1].add_time(time)
        
nodes = [2, 4, 8, 16, 32, 64]
matrices = []
timings = []
size = 0
idx = 0
for n in nodes:
    for f in glob.glob("../runscripts/%s/spmbvs_N%d*.out"%(computer, n)):
        file = open(f)
        
        for line in file:
            if "Testing" in line:
                name = ((line.rsplit("\n")[0]).rsplit("/"))[-1]
                if name in matrices:
                    idx = matrices.index(name)
                else:
                    idx = len(matrices)
                    matrices.append(name)
                    timings.append(Matrix(name))
                timings[idx].add_test(n)
            elif "Time" in line:
                name = (line.rsplit(':'))[0]
                time = (float)(((line.rsplit('\n')[0]).rsplit(' '))[-1])
                timings[idx].add_time(name, time)
        file.close()


print(timings[-1].name)
print(timings[-2].std[-1].time)
print(timings[-2].early_csc[-1].time)


## Plot all matrices
xdata = np.arange(len(matrices))
std = [timings[i].std[-1].time for i in xdata]
csc = [timings[i].csc[-1].time for i in xdata]
early_csc = [timings[i].early_csc[-1].time for i in xdata]
speedup_csc = [std[i] / csc[i] for i in xdata]
speedup_early_csc = [std[i] / early_csc[i] for i in xdata]
# Sort all lists by speedup_early_csc
combined = sorted(zip(matrices, speedup_csc, speedup_early_csc), key=lambda x: x[2])
matrices_sorted, speedup_csc_sorted, speedup_early_csc_sorted = zip(*combined)
print(matrices_sorted, speedup_csc_sorted, speedup_early_csc_sorted )
plt.add_luke_options()
plt.line_plot(speedup_csc_sorted, xdata, color='blue', label = "CSC")
plt.line_plot(speedup_early_csc_sorted, xdata, color='red', label = "Early CSC")
plt.add_labels("Matrix", "Speedup")
plt.add_anchored_legend(ncol=2)
plt.set_xticks(xdata, matrices_sorted, rotation='vertical', fontsize=10)
plt.save_plot("dane_n64_speedup.pdf")
plt.plt.clf()


## Plot Scaling
xdata = np.arange(len(nodes))
std = [timings[-2].std[i].time for i in xdata]
csc = [timings[-2].csc[i].time for i in xdata]
early_csc = [timings[-2].early_csc[i].time for i in xdata]
speedup_csc = [std[i] / csc[i] for i in xdata]
speedup_early_csc = [std[i] / early_csc[i] for i in xdata]
plt.add_luke_options()
plt.line_plot(std, xdata, color='black', label = "Standard")
plt.line_plot(csc, xdata, color='blue', label = "CSC")
plt.line_plot(early_csc, xdata, color='red', label = "Early CSC")
plt.add_labels("Nodes", "Time (seconds)")
plt.set_scale('linear', 'log')
plt.add_anchored_legend(ncol=3)
plt.set_xticks(xdata, nodes)
plt.save_plot("dane_%s.pdf"%matrices[-2][:-3])
plt.plt.clf()


plt.add_luke_options()
plt.line_plot(speedup_csc, xdata, color='blue', label = "CSC")
plt.line_plot(speedup_early_csc, xdata, color='red', label = "Early CSC")
plt.add_labels("Nodes", "Speedup")
plt.add_anchored_legend(ncol=2)
plt.set_xticks(xdata, nodes)
plt.save_plot("dane_%s_speedup.pdf"%matrices[-2][:-3])
plt.plt.clf()
