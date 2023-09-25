from billiard import *
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
#https://towardsdatascience.com/how-to-animate-plots-in-python-2512327c8263



def momentum_plot(m_min=1, m_max=5, contour=False):
    X = np.arange(0.01001, 1, 0.02)
    M = np.arange(m_min + 0.00001, m_max, 0.02)
    Z = [[mom3(x, m3=m) for x in X] for m in M]
    X, M = np.meshgrid(X, M)
    Z = np.array(Z).reshape(X.shape)

    # would like to find a way to split the data into continuous
    # pieces so that we see separate sheets ... but I am not sure
    # how well plot_surface would handle that.
    if contour:
        plt.contour(X,M,Z)
        plt.xlabel('x1')
        plt.ylabel('mass 3')
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.plot_surface(X, M, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.set_xlabel('x1')
        ax.set_ylabel('mass 3')
        ax.set_zlabel('final momentum')
        ax.set_zlim(0.8, 1.01)
    plt.show()


def time_plot(m_min=2.001, m_max=8.001):
    # 3d plot of time to complete (shows the reciprocal, because
    # time -> infinity near break points)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    X = np.arange(0.01001, 1, 0.02)
    M = np.arange(m_min + 0.00001, m_max, 0.02)
    Z = [[1/time3(x, m3=m) for x in X] for m in M]
    # Z = [[min(10, time3(x, m3=m)) for x in X] for m in M]
    X, M = np.meshgrid(X, M)
    Z = np.array(Z).reshape(X.shape)

    ax.plot_surface(X, M, Z,  cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim(0, 0.6)
    plt.show()


def position_plot(x = 0.5, mass2=1, mass3=2, window=None, phase=False):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,6) # does this do anything?
    result = run_bounces([[1, x, 0], [mass2, 1, 0], [mass3, 1.1, -1]])
    paths, time = result['locations'], result['times']
    tmax = time[-1]
    xmax = paths[-1][-1]
    if window is not None:
        tmax, xmax = window
        dt = tmax - time[-1]
        if dt > 0:
            dt1 = time[-1] - time[-2]
            last, previous = paths[-1], paths[-2]
            for i in range(len(last)):
                last[i] += dt * (last[i] - previous[i]) / dt1
            time[-1] = tmax
    _pos_plot(ax, time, paths, x, mass2, mass3, tmax, xmax)
    plt.show()


def _pos_plot(ax, time, paths, x, mass2, mass3, tmax, xmax):
    ax.plot(time, paths)
    ax.set_title("x1 = {0:6.3f}, m2 = {1:6.3f}, m3 = {2:6.3f}".format(x, mass2, mass3))
    ax.set_xlabel('time')
    ax.set_ylabel('position')
    ax.set_xlim([0, tmax])
    ax.set_ylim([0, xmax])


def position_animation(mass2=1, mass3=2, xmin=0.0001, xmax=0.9999, save=False, window=None, delay=200, phase=False):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,6)
    xmin = min(max(0.0001, xmin), 0.95)
    xmax = min(0.9999, max(xmax, xmin+0.05))
    N = math.ceil(100 / math.sqrt(1/(xmax-xmin)))
    dx = (xmax-xmin)/N
    print(xmin, xmax, N, dx)
    X = np.arange(xmin, xmax+dx/2, dx)
    plot_data, vel_data = [], []
    tmax = 0
    xmax = 0
    for x in X:
        result = run_bounces([[1, x, 0], [mass2, 1, 0], [mass3, 1.1, -1]])
        path, time = result['locations'], result['times']
        if phase:
            vel = result['vel']
            vel_data.append([ [v[0] for v in vel], [v[1] for v in vel] ])
        tmax = max(tmax, time[-1])
        xmax = max(xmax, path[-1][-1])
        plot_data.append((path, time))
    if window is not None:
        tmax, xmax = window
    # extend lines so that all plots fill up the whole time
    for path, time in plot_data:
        dt = tmax - time[-1]
        if dt == 0: continue
        dt1 = time[-1] - time[-2]
        last, previous = path[-1], path[-2]
        for i in range(len(last)):
            last[i] += dt * (last[i] - previous[i]) / dt1
        time[-1] = tmax
    if phase:
        ax2 = fig.add_axes([0.2, 0.5, 0.3, 0.3])
        ax2.axis('off')

    def animate_func(num):
        ax.clear()
        _pos_plot(ax, plot_data[num][1], plot_data[num][0], X[num], mass2, mass3, tmax, xmax)
        if phase:
            ax2.clear()
            phase_plane(X[num], mass2, mass3, ax2, vel_data[num])

    ani = animation.FuncAnimation(fig, animate_func,
                                interval=delay, frames=len(plot_data),
                                repeat=False, repeat_delay=500)
    plt.show()
    if save:
        print("saving animation ...")
        FFwriter = animation.FFMpegWriter(fps=10)
        ani.save('animation.mp4', writer = FFwriter)
        f = f"position plot-m2={mass2}, m3={mass3}.gif"
        writergif = animation.PillowWriter(fps=10)
        ani.save(f, writer=writergif)


def two_ball_animation(save=False, window=None):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,6)

    M = [ x/25 for x in range(26,200) ] \
        + [ x/10 for x in range(80,120) ] \
        + [ x/5 for x in range(61,100) ] \
        + [ x/2 for x in range(41,80) ] \
        + [ x for x in range(41,100) ]
        
    plot_data = []
    tmax = 0
    xmax = 0
    for m in M:
        result = run_bounces([[1, 1, 0], [m, 2, -1]])
        path, time = result['locations'], result['times']
        tmax = max(tmax, time[-1])
        xmax = max(xmax, path[-1][-1])
        plot_data.append((path, time))
    if window is not None:
        tmax, xmax = window
    # extend lines so that all plots fill up the whole time
    for path, time in plot_data:
        dt = tmax - time[-1]
        if dt == 0: continue
        dt1 = time[-1] - time[-2]
        last, previous = path[-1], path[-2]
        for i in range(len(last)):
            last[i] += dt * (last[i] - previous[i]) / dt1
        time[-1] = tmax
   
    def animate_func(num):
        ax.clear()
        ax.plot(plot_data[num][1], plot_data[num][0])
        ax.set_title("larger mass = {:6.2f}".format(M[num]))
        ax.set_xlabel('time')
        ax.set_ylabel('position')
        ax.set_xlim([0, tmax])
        ax.set_ylim([0, xmax])

    ani = animation.FuncAnimation(fig,
                                animate_func,
                                interval=50,   
                                frames=len(M),
                                repeat=False,
                                repeat_delay=100)
    plt.show()
    if save:
        print("saving animation ...")
        f = "two-ball animation.gif"
        _writer = animation.PillowWriter(fps=18)
        ani.save(f, writer=_writer)


def all_phases(mass2=1, mass3=2.01):
    breaks = find_breaks(mass2, mass3) + [1]
    X = [[breaks[0]/2]]
    for i in range(len(breaks)-1):
        b1, b2 = breaks[i:i+2]
        X.append([b1*0.999+b2*0.001, b1*0.001+b2*0.999])
    fig, axs = plt.subplots(2,len(X))
    fig.suptitle(f"mass 2 = {mass2:6.3f}, mass 3 = {mass3:6.3f}")
    for i, col in enumerate(X):
        for j, x in enumerate(col):
           axs[j][i].set_title(f"x = {x:5.3f}")
           phase_plane(x, mass2, mass3, axs[j][i])
    axs[1][0].axis('off')
    plt.show()
    # phase_plane(X, mass2, mass3)


def phase_plane(X, mass2=1, mass3=2, ax=None, data=None):
    if not isinstance(X, list):
        X = [X]
    endV1, endV2 = [], []
    show = False
    if ax is None:
        fig, ax = plt.subplots()
        show = True
    # fig.set_size_inches(10,8)
    ax.set_aspect(1)
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])

    ax.add_artist(plt.Circle((0,0), 1, fill=False))
    if data is None:  
        for x in X:
            vel = run_bounces([[1, x, 0], [mass2, 1, 0], [mass3, 1.1, -1]])['vel']
            v1, v2 = [v[0] for v in vel], [v[1] for v in vel]
            ax.plot(v1, v2)
            endV1.append(v1[-1])
            endV2.append(v2[-1])
    else:
        v1, v2 = data
        ax.plot(v1, v2)
        endV1, endV2 = [v1[-1]], [v2[-1]]
    ax.scatter(endV1, endV2)
    ax.axis('off')
    if show: plt.show()

def phase_movie(mass2=1, mass3=2, xmin=0.001, xmax=0.999):
    data = []
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])
    dx = (xmax-xmin)/50
    X = np.arange(xmin, xmax+dx/2, dx)

    for x in X:
        result = run_bounces([[1, x, 0], [mass2, 1, 0], [mass3, 1.1, -1]])
        vel = result['vel']
        v1, v2 = [v[0] for v in vel], [v[1] for v in vel]
        data.append([x, v1, v2])
    def animate_func(num):
        ax.clear()
        ax.set_aspect(1)
        x, v1, v2 = data[num]
        ax.set_xlim([-1.01, 1.01])
        ax.set_ylim([-1.01, 1.01])
        ax.add_artist(plt.Circle((0,0), 1, fill=False))
        ax.axis('off')
        ax.plot(v1, v2)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        plt.scatter([v1[-1]], [v2[-1]])
        ax.set_title("x = {0:5.3f}, m2 = {1:6.3f}, m3 = {2:6.3f}".format(x, mass2, mass3))
 
    ani = animation.FuncAnimation(fig,
                                animate_func,
                                interval=250,   
                                frames=len(data),
                                repeat=False,
                                repeat_delay=800)
    plt.axis('off')
    plt.show()


def show_bounces(x = 0.5, mass2=1, mass3=2.00001, save = False):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,6)

    bounces = run_bounces([[1, x, 0], [mass2, 1, 0], [mass3, 1.1, -1]])['bounces']
    X = bounces[0]
    index, t, dt, time = 1, 0.02, 0.02, 0.0
    data = [[time, X]]
    vel = bounces[1][1]
    tmax = bounces[1][0]
    while True:
        if t < tmax:
            X = [X[i] + vel[i] * dt for i in range(3)]
        else:
            X = [X[i] + vel[i] * (dt - (t - tmax)) for i in range(3)]
            index += 1
            if index == len(bounces): break
            vel = bounces[index][1]
            t -= tmax
            X = [X[i] + vel[i] * t for i in range(3)]
            tmax = bounces[index][0]
        if X[-1] > 2:
            dt = X[-1]/100
        time += dt
        t += dt
        data.append([time, X])

    def animate_func(num):
        ax.clear()
        t, [x1, x2, x3] = data[num]
        ax.scatter([x1], [0.1])
        ax.scatter([x2], [0.1])
        ax.scatter([x3], [0.1])
        xmax = max(x3*1.05, 1.3)
        ax.set_xlim([0, xmax])
        ax.set_ylim([0, 0.2])
        ax.text(xmax*0.45, 0.15, "t = {0:6.3f}".format(t))
        ax.set_yticks([])
        ax.spines[['right', 'top']].set_visible(False)

 
    ani = animation.FuncAnimation(fig,
                                animate_func,
                                interval=50,   
                                frames=len(data),
                                repeat=False,
                                repeat_delay=800)
    plt.show()
    if save:
        print("saving animation ...")
        f = "bounce movie.gif"
        _writer = animation.PillowWriter(fps=18)
        ani.save(f, writer=_writer)
            
# show_bounces(mass3=6.01, save=False) # saved
# position_plot(0.5, 1, 6.01) # saved

# show_bounces(x = 0.9, mass2=5, mass3=10)
# position_plot(0.9, 5, 10)
# position_animation(5, 10, window = (5,4))

# position_plot(x=0.1, mass2=1, mass3=4.1)
# position_plot(x=0.5, mass2=1, mass3=4.1, window = (5,4))
# position_animation(1, 4.1, window = (5,4))
# position_animation(1, 9.1, window = (5,4))

# these four plots show the impact of a very small change in
# the initial position of mass1.
# position_plot(x=0.284, mass2=5, mass3=10, window=(5, 4))
# position_plot(x=0.285, mass2=5, mass3=10, window=(5, 4))
# position_plot(x=0.455, mass2=5, mass3=10, window=(5, 4))
# position_plot(x=0.456, mass2=5, mass3=10, window=(5, 4))

# momentum_plot(2, 10, contour=False)
# position_animation(mass2=5, mass3=10, window=(5, 4), save=True)
# position_animation(mass2=1.3, mass3=3, window=(5, 4), save=True)
# position_animation(mass3=5.828, window=(5, 4), save=True)
position_animation(mass3=6.5, window=(4, 3), save=True, phase=True, delay=100)

# two_ball_animation(window = (8, 6), save=False)

# phase_plane(x=0.1, mass3=4.1)
# phase_plane(x=0.22, mass3=4.1)
# phase_plane([0.22, 0.5, 0.75], mass3=4.1)
# phase_plane(x=0.75, mass3=4.1)
# position_plot(0.75, 1, 4.1)
# position_animation(1, 3.1, window=(5,4))
# all_phases(mass3=5.1)
# all_phases(mass3=8.1)
# all_phases(mass3=10.1)
# print(find_breaks(mass2=2, mass3=8.1))
# momentum_plot()
# all_phases(mass2=2, mass3=8.1)
# 0.04007999420266015, 0.15026278495888573  yes
# 0.42864721298317765, 0.8278688430796131   yes
# 0.923217544556664, 0.9505882263193594]

# position_animation(mass2=2, mass3=8.1, xmin=0.0401, xmax=0.1502, window=(3,2), delay=500)
# position_animation(mass2=2, mass3=8.1, xmin=0.4287, xmax=0.8278, window=(3,2), delay=500)
# phase_plane([0.1, 0.5, 0.9], mass3=6.1)
# print(find_breaks(mass3=10.0001))
# phase_movie(mass3=10.0001)
# phase_movie(mass3=10.0001, xmin=0, xmax=0.301)
# phase_movie(mass3=10.0001, xmin=0.302, xmax=0.379)
# phase_movie(mass3=10.0001, xmin=0.380, xmax=0.712)
# phase_movie(mass3=10.0001, xmin=0.713, xmax=0.949)
# phase_movie(mass3=10.0001, xmin=0.950, xmax=1)
# position_animation(mass3=10.0001, window=(3,2), xmin=0.95)
# print(find_breaks(mass3=5.0001))
# phase_movie(mass3=5.0001)
# phase_movie(mass3=5.0001, xmin=0, xmax=0.428)
# phase_movie(mass3=5.0001, xmin=0.429, xmax=0.891)
# phase_movie(mass3=5.0001, xmin=0.892, xmax=1)


# Hypothesis: It looks like, given the "breaks" b1, b2, ...
# for a particular configuration:
#  * there is only one (simple) phase plane for x<b1.
#  * between 
# phase_movie(mass3=8.001)
# print(run_bounces([[1, 0.5, 0], [1, 1, 0], [2, 1.1, -1]])['bounces'])
# show_bounces()

# momentum_plot(2,10)
# time_plot(2, 10)
# print([mom3(x/10, 1, 5.3) for x in range(1,10)])
# X = [x/10 for x in range(1,10)] # np.arange(0.01001, 1, 0.1)
# M = [5.3, 8.3]
# Z = [[mom3(x, 1, m) for x in X] for m in M]
# print(Z)


# X = np.arange(1, 5, 1)
# M = np.arange(-2, 2, 1)
# Z = [[x+m/5 for x in X] for m in M]
# X, M = np.meshgrid(X, M)
# print(X)
# print(M)
# print(Z)
# X1 = np.reshape(X, (2,8))
# print(X1)
# # X[0] = X[0][1:3]
# # Z[0