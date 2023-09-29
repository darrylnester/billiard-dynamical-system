from billiard import *
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm
import matplotlib as mpl

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "b", "k", "g", "#380282"]) 

# reference: https://towardsdatascience.com/how-to-animate-plots-in-python-2512327c8263


def mass_text(masses):
    return 'masses ' + make_title([m[0] for m in masses])
    # txt = 'masses ' + ', '.join([f"{m[0]:6.3f}" for m in masses])
    # return txt.replace('.000', '')


def three_decimals(x, no_zeros=True, no_trail=False):
    txt = f"{x:10.3f}".strip()
    if no_trail:
        while txt[-1] == '0':
            txt = txt[0:-1]
        if txt[-1] == '.':
            txt = txt[0:-1]
    if no_zeros:
        txt = txt.replace('.000', '')
    return txt


def make_title(nums, names=None, no_zeros=True):
    if names is None:
        txt = ', '.join([three_decimals(n, no_trail=True) for n in nums])
    else:
        txt = ', '.join([names[i]+'='+three_decimals(nums[i], no_trail=True) for i in range(len(nums))])
    return "".join(txt)
#    return " ".join(txt.replace('.000', '').strip().split())


def momentum_plot(m_min=1.5, m_max=5, contour=False, azim_range=None, elev_range=None, delay=100, save=False):
    X = np.arange(0.01001, 1, 0.01)
    M = np.arange(m_min + 0.00001, m_max, 0.02)
    Z = [[mom3(x, m3=m) for x in X] for m in M]
    X, M = np.meshgrid(X, M)
    Z = np.array(Z).reshape(X.shape)
    # would like to find a way to split the data into continuous
    # pieces so that we see separate sheets ... but I am not sure
    # how well plot_surface would handle that.
    if contour:
        plt.contour(X,M,Z)
        plt.xlabel('x')
        plt.ylabel('m3')
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        fig.set_size_inches(6,5)
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        ax.plot_surface(X, M, Z,
                        cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.set_xlabel('x')
        ax.set_ylabel('m3')
        ax.set_zlabel('final momentum')
        ax.set_zlim(0.8, 1.01)
        if azim_range is not None or elev_range is not None:
            if azim_range is None:
                azim_range = [45] * len(elev_range)
            if elev_range is None:
                elev_range = [30] * len(azim_range)
            def animate_func(num):
                ax.view_init(elev=elev_range[num], azim=azim_range[num])
            ani = animation.FuncAnimation(fig, animate_func,
                                interval=delay, frames=len(azim_range),
                                repeat=False, repeat_delay=500)
    plt.show()
    save_check(ani, save, f"momentum plot-{m_min} to {m_max}", 10)


def save_check(ani, save, title="animation", frames_per_sec=10):
    if save is not False:
        print("saving animation ...")
        if save == "mp4":
            title += ".mp4"
            ani_writer = animation.FFMpegWriter(fps=frames_per_sec)
        else:
            title += ".gif"
            ani_writer = animation.PillowWriter(fps=frames_per_sec)
        ani.save(title, writer = ani_writer)


def time_plot(m_min=2.001, m_max=8.001):
    # 3d plot of time to complete (shows the reciprocal, because
    # time -> infinity near break points)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    fig.set_size_inches(8,6)
    X = np.arange(0.01001, 1, 0.02)
    M = np.arange(m_min + 0.00001, m_max, 0.02)
    Z = [[1/time3(x, m3=m) for x in X] for m in M]
    X, M = np.meshgrid(X, M)
    Z = np.array(Z).reshape(X.shape)

    ax.plot_surface(X, M, Z,  cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim(0, 0.6)
    plt.show()


def _pos_plot(ax, time, paths, mass_list, tmax, xmax, ax2, vel_data):
    ax.plot(time, paths)
    ax.set_title(mass_text(mass_list))
    ax.set_xlabel('time')
    ax.set_ylabel('position')
    ax.set_xlim([0, tmax])
    ax.set_ylim([0, xmax])
    if ax2 is not None:
        ax2.clear()
        phase_plane(masses=mass_list, ax=ax2, data=vel_data)


def _norm_vel(bounce_results):
    vel, c = bounce_results['vel'], bounce_results['coeffs']
    return [ [v[0]*c[0] for v in vel], [v[1]*c[1] for v in vel] ] 


def position_plot(x = 0.5, m2=1, m3=2, masses=None, window=None, phase=False, phase_location=[0.1, 0.5, 0.35, 0.35]):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,6)
    if masses is None:
        masses = [[1, x, 0], [m2, 1, 0], [m3, 1.1, -1]]
    results = run_bounces(masses)
    paths, time = results['locations'], results['times']
    tmax = time[-1]
    xmax = paths[-1][-1]
    ax2, vel_data = None, None
    if phase:
        vel_data = _norm_vel(results)
        ax2 = fig.add_axes(phase_location)
        ax2.axis('off')

    if window is not None:
        tmax, xmax = window
        dt = tmax - time[-1]
        if dt > 0:
            dt1 = time[-1] - time[-2]
            last, previous = paths[-1], paths[-2]
            for i in range(len(last)):
                last[i] += dt * (last[i] - previous[i]) / dt1
            time[-1] = tmax
    _pos_plot(ax, time, paths, masses, tmax, xmax, ax2, vel_data)
    plt.show()


def position_animation(m2=1, m3=2, xmin=0.0001, xmax=0.9999, masses=None, save=False, window=None, delay=200, phase=False, phase_location=[0.1, 0.5, 0.35, 0.35]):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,6)
    xmin = min(max(0.0001, xmin), 0.95)
    xmax = min(0.9999, max(xmax, xmin+0.05))
    N = math.ceil(100 / math.sqrt(1/(xmax-xmin)))
    dx = (xmax-xmin)/N
    X = np.arange(xmin, xmax+dx/2, dx)
    plot_data, vel_data = [], []
    tmax, xmax = 0, 0
    if masses is None: 
        masses = [[1, 0.1, 0], [m2, 1, 0], [m3, 1.1, -1]]
    for x in X:
        masses[0][1] = x
        results = run_bounces(masses)
        path, time = results['locations'], results['times']
        tmax = max(tmax, time[-1])
        xmax = max(xmax, path[-1][-1])
        plot_data.append( (path, time) )
        vel_data.append( _norm_vel(results) )
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
    ax2 = None
    if phase:
        ax2 = fig.add_axes(phase_location)
        ax2.axis('off')

    def animate_func(num):
        ax.clear()
        _pos_plot(ax, plot_data[num][1], plot_data[num][0], masses, tmax, xmax, ax2, vel_data[num])

    ani = animation.FuncAnimation(fig, animate_func,
                                interval=delay, frames=len(plot_data),
                                repeat=False, repeat_delay=500)
    plt.show()
    save_check(ani, save, "position animation-" + mass_text(masses), 10)


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
    save_check(ani, save, "two-ball animation", 18)


def all_phases(m2=1, m3=2.01):
    breaks = find_breaks(m2, m3) + [1]
    X = [[breaks[0]/2]]
    for i in range(len(breaks)-1):
        b1, b2 = breaks[i:i+2]
        X.append([b1*0.999+b2*0.001, b1*0.001+b2*0.999])
    fig, axs = plt.subplots(2,len(X))
    fig.suptitle(make_title([m2, m3], ['mass 2', 'mass 3']))
        # f"mass 2 = {m2:6.3f}, mass 3 = {m3:6.3f}")
    for i, col in enumerate(X):
        for j, x in enumerate(col):
           axs[j][i].set_title(make_title([x],['x']))
           phase_plane(x, m2, m3, ax=axs[j][i])
    axs[1][0].axis('off')
    plt.show()
 

def phase_plane(X=0.1, m2=1, m3=2, masses=None, ax=None, data=None):
    # creates one or more phase plane plots. If data is
    # given, use it; otherwise, find the velocities by
    # running the bounces for the given mass configuration.
    if not isinstance(X, list):
        X = [X]
    endV1, endV2 = [], []
    show = False
    if ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(5,5)
        plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        show = True
    ax.axis('off')
    ax.set_aspect(1)
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])

    ax.add_artist(plt.Circle((0,0), 1, fill=False))
    if data is None:  
        if masses is None:
            masses = [[1, 0.1, 0], [m2, 1, 0], [m3, 1.1, -1]]
        for x in X:
            masses[0][1] = x
            results = run_bounces(masses)
            v1, v2 = _norm_vel(results)
            ax.plot(v1, v2, ':', linewidth=2)
            endV1.append(v1[-1])
            endV2.append(v2[-1])
    else:
        v1, v2 = data
        ax.plot(v1, v2, ':', linewidth=2)
        endV1, endV2 = [v1[-1]], [v2[-1]]
    ax.scatter(endV1, endV2)
    if show: plt.show()


def phase_movie(m2=1, m3=2, xmin=0.001, xmax=0.999):
    data = []
    fig, ax = plt.subplots()
    fig.set_size_inches(5,5)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    ax.set_aspect(1)
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])
    dx = (xmax-xmin)/50
    X = np.arange(xmin, xmax+dx/2, dx)

    for x in X:
        result = run_bounces([[1, x, 0], [m2, 1, 0], [m3, 1.1, -1]])
        v1, v2 = _norm_vel(result)
        data.append([x, v1, v2])
    def animate_func(num):
        ax.clear()
        ax.set_aspect(1)
        x, v1, v2 = data[num]
        ax.set_xlim([-1.01, 1.01])
        ax.set_ylim([-1.01, 1.01])
        ax.add_artist(plt.Circle((0,0), 1, fill=False))
        ax.axis('off')
        ax.plot(v1, v2, ':', linewidth=2)
        plt.scatter([v1[-1]], [v2[-1]])
        ax.set_title(make_title([x,m2,m3],['x','m2','m3']))
        # ax.set_title(f"x = {x:5.3f}, m2 = {m2:6.3f}, m3 = {m3:6.3f}")
 
    ani = animation.FuncAnimation(fig,
                                animate_func,
                                interval=250,   
                                frames=len(data),
                                repeat=False,
                                repeat_delay=800)
    plt.axis('off')
    plt.show()


def show_bounces(x = 0.5, m2=1, m3=2.00001, masses=None, min_x=1.3, save = False, show_vel=True):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(8,3)
    plt.subplots_adjust(left=0.03, right=0.97, top=0.99, bottom=0.08)

    if masses is None:
        masses = [[1, x, 0], [m2, 1, 0], [m3, 1.1, -1]]
        # note: the last mass must be the furthest from the wall

    title = mass_text(masses)
    N = len(masses)
    bounces = run_bounces(masses)['bounces']    
    X = bounces[0]
    index, time = 1, 0.0
    data = [[time, X, 0]]
    tmax, vel = bounces[1]
    dt = min(0.005, tmax/2)
    t = dt
    while True:
        if t < tmax:
            X = [X[i] + vel[i] * dt for i in range(N)]
        else:
            X = [X[i] + vel[i] * (dt - (t - tmax)) for i in range(N)]
            index += 1
            if index == len(bounces): break
            vel = bounces[index][1]
            t -= tmax
            X = [X[i] + vel[i] * t for i in range(N)]
            tmax = bounces[index][0]
            dt = min(0.005, tmax/2.5)
        if X[-1] > 2:
            dt = X[-1]/100
        time += dt
        t += dt
        data.append([time, X, index-1])

    def animate_func(num):
        ax.clear()
        t, X, c = data[num]
        for x in X:
            ax.scatter([x], [0.13])
        xmax = max(X[-1]*1.05, min_x)
        ax.set_xlim([0, xmax])
        ax.set_ylim([0, 0.25])
        ax.text(xmax/2, 0.23, "t = {0:6.3f}".format(t), horizontalalignment='center')
        ax.text(xmax/2, 0.19, title, horizontalalignment='center')
        ax.text(xmax/2, 0.16, f"{c} collisions", horizontalalignment='center')
        if show_vel or num==len(data)-1:
            vel = bounces[c+1][1] # bounces[-1][1]
            last_position = -xmax/50
            for i, x in enumerate(X):
                pos = max(x-xmax/160, last_position + xmax/50)
                ax.text(pos, 0.12, f"v{i+1} = {vel[i]:8.5f}", rotation=90, verticalalignment='top')
                last_position = pos
        ax.set_yticks([])
        ax.spines[['right', 'top']].set_visible(False)
 
    ani = animation.FuncAnimation(fig,
                                animate_func,
                                interval=20,   
                                frames=len(data),
                                repeat=False,
                                repeat_delay=800)
    plt.show()
    if save is not False:
        print("saving animation ...")
        if save == "mp4":
            f = f"bounce movie-{title}.mp4"
            ani_writer = animation.FFMpegWriter(fps=18)
        else:
            f = f"bounce movie-{title}.gif"
            ani_writer = animation.PillowWriter(fps=18)
        ani.save(f, writer = ani_writer)


def phase_mass_movie(x=0.5, m3=(1.0001,5), save=False, frames=100):
    m3_min, m3_max = m3
    data = []
    fig, ax = plt.subplots()
    fig.set_size_inches(5.5,5)
    plt.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.01)
    ax.set_aspect(1)
    ax.set_xlim([-1.01, 1.01])
    ax.set_ylim([-1.01, 1.01])
    for m in np.arange(m3_min, m3_max, (m3_max-m3_min)/frames):
        result = run_bounces([[1, x, 0], [1, 1, 0], [m, 1.1, -1]])
        v1, v2 = _norm_vel(result)
        data.append([m, v1, v2])
    def animate_func(num):
        ax.clear()
        ax.set_aspect(1)
        m, v1, v2 = data[num]
        ax.set_xlim([-1.01, 1.01])
        ax.set_ylim([-1.01, 1.01])
        ax.add_artist(plt.Circle((0,0), 1, fill=False))
        ax.axis('off')
        ax.plot(v1, v2, ':', linewidth=2)
        plt.scatter([v1[-1]], [v2[-1]])
        ax.set_title(f"x={x}, m3={three_decimals(m,False)}")
 
    ani = animation.FuncAnimation(fig,
                                animate_func,
                                interval=250,   
                                frames=len(data),
                                repeat=False,
                                repeat_delay=800)
    plt.axis('off')
    plt.show()
    save_check(ani, save, f"phase-mass-movie-{three_decimals(m3_min,no_trail=True)} to {three_decimals(m3_max,no_trail=True)}", 10)


###### CALLS TO THE ABOVE FUNCTIONS TO PRODUCE PLOTS ######


# show_bounces(m3=6.01) # saved
# show_bounces(masses = [[1, 0.5, 0], [20, 1, -1]])
# show_bounces(masses = [[1, 0.1, 0], [2, 0.3, 0], [4, 0.5, 0], [8, 1, -1]])
# show_bounces(masses = [[1, 0.1, 0], [2, 0.2, 0], [3, 0.3, 0], [4, 0.4, 0], [5, 0.5, 0], [6, 0.6, -1]], save='mp4')
# show_bounces(masses = [[1, 0.1, 0], [1, 0.2, 0], [2, 0.3, 0], [3, 0.4, -1]], save='mp4')
# position_plot(masses = [[1, 0.1, 0], [3, 0.2, 0], [2, 0.3, 0], [3, 0.4, -1]])
# show_bounces(masses = [[1, 0.5, 0], [100, 0.55, -1]], min_x=0.6, show_vel=False, save="mp4")
# show_bounces(masses = [[1, 0.4, 0], [10000, 0.5, -1]], min_x=0.6, show_vel=False, save="mp4")

# position_plot(0.5, 1, 6.01) # saved


# position_animation(m2=2, m3=7.5, window = (5,4), phase=True)
# position_animation(m3=5.828, window = (5,4), phase=True, save='mp4')
# position_animation(m2=1.3, m3=3, window = (4,2.5), phase=True, save='mp4')

# position_animation(m2=5, m3=10, window = (5,3.5), phase=True, save='mp4')
# position_animation(3, 5, window = (5,4), phase=True)

# these four plots show the impact of a very small change in
# the initial position of mass1.
# position_plot(x=0.284, m2=5, m3=10, window=(5, 4), phase=True)
# position_plot(x=0.285, m2=5, m3=10, window=(5, 4), phase=True)
# position_plot(x=0.455, m2=5, m3=10, window=(5, 4), phase=True)
# position_plot(x=0.456, m2=5, m3=10, window=(5, 4), phase=True)

# position_animation(m2=5, m3=10, window=(5, 4), phase=True, save='mp4')
# position_animation(m2=1.3, m3=3, window=(5, 4), phase=True, save='mp4')
# position_animation(m3=5.828, window=(5, 4), phase=True, save='mp4')
# position_animation(m2=1, m3=6, window=(4.5, 3), save='mp4', phase=True)

# two_ball_animation(window = (8, 6), save=False)

# phase_plane([0.1, 0.5, 0.9], m3=6.1)
# phase_movie(1, 5)
# phase_plane(0.7, 1, 5)
# position_plot(0.7, 1, 5, window=(4,3), phase=True, phase_location=[0.1, 0.4, 0.45, 0.45])
# phase_mass_movie(x=0.7, m3=(1.0001,8), save='mp4', frames=210)
# M = np.arange(1.0001,8, 0.005)
# L = [mom3(0.6, 1, m) for m in M]
# plt.plot(M,L)
# plt.show()

# phase_plane([0.1, 0.5], m3=4.1)
# position_plot(0.75, 1, 4.1)
# position_animation(1, 3.1, window=(5,4))
# all_phases(m3=5.1)
# all_phases(m3=8.1)
# all_phases(m3=10.1)
# print(find_breaks(m2=2, m3=8.1))

# position_animation(m2=2, m3=8.1, xmin=0.0401, xmax=0.1502, window=(3,2), delay=500)
# position_animation(m2=2, m3=8.1, xmin=0.4287, xmax=0.8278, window=(3,2), delay=500)

# azim = np.arange(-80,0,0.8)
# elev = np.arange(20,90,70/len(azim))
# # momentum_plot(1.5,10)
# momentum_plot(1.5,10, azim_range=azim, elev_range=elev, delay=120, save='mp4')
# momentum_plot()
# time_plot(2, 10)

