import glob
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
from time import sleep
from matplotlib.pyplot import show, plot
import sys

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

if (len(sys.argv)!=2):
    print "Wrong arg number:"+str(len(sys.argv))+" for script, python plot.py ./your_file_path"
    sys.exit(-1)
print sys.argv
fig = plt.figure()
files = glob.glob(sys.argv[1])
files.sort(key= natural_keys)

for i in xrange(0,  len(files)):
    file_name = files[i]
    flag = 0
    with open(file_name) as f:
        for line in f:
            numbers_float = map(float, line.split())
            if len(numbers_float)==3:
                if flag == 0:
                    ax = fig.add_subplot(111, projection='3d')
                    ax.set_xlabel('X')
                    ax.set_ylabel('Y')
                    ax.set_zlabel('Z')
                    ax.view_init(50, 40)
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    ax.set_zlim(0, 1)
                    ax.set_title(file_name+"("+str (i+1)+"/"+str(len(files))+")")
                    flag = 1
                ax.scatter(numbers_float[0], numbers_float[1], numbers_float[2])
            if len(numbers_float) ==2:
                if flag == 0:
                    ax = fig.add_subplot(111)
                    ax.set_xlabel('X')
                    ax.set_ylabel('Y')
                    ax.set_xlim(0,1)
                    ax.set_ylim(0,1)
                    ax.set_title(file_name+"("+str (i+1)+"/"+str(len(files))+")")
                    flag = 1
                ax.scatter(numbers_float[0], numbers_float[1])
        f.close()
    plt.draw()

    if i != len(files)-1 :
        plt.show(block=False)   #this creates an empty frozen window.
        sleep(0.01)
        ax.cla()
    else:
        plt.show(block=True)

