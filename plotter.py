import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

data = np.genfromtxt('testData.dat', delimiter=',')


fig = plt.figure()
ax = plt.axes(xlim=(-1,1), ylim=(0,2))
line, = ax.plot([], [], lw=2)

def init ():
    line.set_data([],[])
    return line,

def animate(i):
    x = np.linspace(-1.0, 1.0,1000)
    y = data[i]
    line.set_data(x,y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1000, interval=20, blit=True)

anim.save('shock_full.mp4', fps=50)

plt.show()
