import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

T = 100
L = 10
# Note, (c*dt/dx) must be less than or equal to 1, otherwise numerical disperation relation blows up :(
c = 2
dt = 0.05
x_space, dx = np.linspace(0, L, num=100, retstep=True)


def spatial_derivative(u, dx):
    out = [0] * (len(u))

    for index in range(1,len(u)):
        out[index] = (u[index] - u[index-1])/dx

    # extrapolate derivative to beginning
    # this isn't a very good approximation but oh well
    # otherwise, we will lose a point each iteration and yeah thats not good.
    out[0] = out[1]

    return np.array(out)

def spatial_second(u, dx):
    out = np.zeros(len(u))
    for i in range(1, len(u)-1):
    
        out[i] = u[i-1] - 2*u[i] + u[i+1]

    return out / (dx*dx)

def time_derivative(u, u_old, dt):
    out = [ u[i] - u_old[i] for  i in range(len(u))]
    return np.array(out)/dt

def main():

    # using flattened normal distribution
    def f(x_space, t):
        out = np.exp(-(x_space-5*L/8.0-c*t)**2)
        return out 
    def g(x_space, t):
        return f(x_space, t)

    u_arr = [
            f(x_space, 0) + g(x_space, 0),
            f(x_space, dt) + g(x_space, -dt),
            f(x_space, 2*dt) + g(x_space, -2*dt)
    ]

    rsq=(c*dt/dx)**2
    t = 0
    i = 2
    while t < T:
        u = u_arr[i]
        prev_u = u_arr[i-1]

        """
        next_u = np.zeros(len(u))
        for a in range(1, len(u)-1): 
            next_u[a]  = 2*(1-rsq)*u[a]-prev_u[a]+rsq*(u[a-1]+u[a+1])
        """

        # compute u_xx
        u_xx = np.zeros(len(u))
        for j in range(1, len(u_xx) -1):
            u_xx[j] = u[j-1] -2*u[j] + u[j+1]

        # wave equation solved for u_t
        next_u = rsq* u_xx  + 2*u - prev_u

        # boundary conditions
        next_u[0] = 0
        next_u[len(u)-1]=0

        u_arr.append(next_u)
        print(next_u)

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
