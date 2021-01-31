import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

T = 10
dt = 0.1
L = 20
x_space, dx = np.linspace(0, L, num=500, retstep=True)
c = 0.4 # speed of wave


def spatial_derivative(u, dx):
    out = [0] * (len(u))

    for index in range(1,len(u)):
        out[index] = (u[index] - u[index-1])/dx

    # extrapolate derivative to beginning
    # this isn't a very good approximation but oh well
    # otherwise, we will lose a point each iteration and yeah thats not good.
    out[0] = out[1]

    return np.array(out)

def time_derivative(u, u_old, dt):
    out = [ u[i] - u_old[i] for  i in range(len(u))]
    return np.array(out)/dt

def main():
    # let u_arr just be sin(x+ct) + sin(x-ct)
    def f(x_space, t):
        return np.exp(-(x_space-L/2-t)**2)
    def g(x_space, t):
        return f(x_space, t)
    u_arr = [
            f(x_space, 0) + g(x_space, 0),
            f(x_space, dt) + g(x_space, -dt),
            f(x_space, 2*dt) + g(x_space, -2*dt)
    ]

    t = 0
    i = 2
    while t < T:
        u = u_arr[i]
        prev_u = u_arr[i-1]
        prev_prev_u = u_arr[i-2]

        # compute derivatives
        u_x = spatial_derivative(u, dx)
        u_xx = spatial_derivative(u_x, dx)

        # finally apply wave equation
        u_tt = (c ** 2) * u_xx

        # now that we know u_tt, find next u
        next_u = (dt ** 2) * u_tt  + 2*prev_u - prev_prev_u
        u_arr.append(next_u)

        """
        print(i)
        print("u_x", u_x)
        print("u_xx", u_xx)
        print("next_u", next_u)
        print("u_prev", u)
        print()
        """
        i += 1
        t += dt

    animate(u_arr)

def animate(u_arr):
    fig, ax = plt.subplots()
    ln, = plt.plot([], [])

    def init():
        ax.set_xlim(0, L)
        ax.set_ylim(-2, 2)
        return ln,

    def update(frame):
        ln.set_data(x_space, u_arr[frame])
        return ln,

    ani = FuncAnimation(fig, update, frames=range(int(T/dt)),
                    init_func=init, blit=True, interval=100, repeat=False)
    plt.show()


def testing():
    t = np.arange(0, 2*np.pi, 0.1)
    plt.plot(t, np.sin(t))
    plt.plot(t, spatial_derivative(np.sin(t), 0.1))
    plt.show()

#testing()
if __name__ == "__main__":
   main()
