import matplotlib.pyplot as plt

dx = lambda alf1, alf3, lamb, m, x, y, z: alf1*z*y - alf3*z*x - lamb*x + m*y
dy = lambda alf0, r0, y, z: r0*y - alf0*z*y
dz = lambda ganm, b, x, y, z: ganm*z*(x + y) - b*z

def RungeKutta(dx, dy, dz, alf1, alf2, alf3, b, ganm, lamb, m, r, x, y, z, h, t, tn, rnd):
    r0 = r - m
    alf0 = alf1 + alf2
    ts = [t]
    xs = [x]
    ys = [y]
    zs = [z]
    while t < tn:
        k1 = dx(alf1, alf3, lamb, m, x, y, z)
        k2 = dx(alf1, alf3, lamb, m, x + h*k1/2, round(y + (h/2), rnd + 1), round(z + (h/2), rnd + 1))
        k3 = dx(alf1, alf3, lamb, m, x + h*k2/2, round(y + (h/2), rnd + 1), round(z + (h/2), rnd + 1))
        k4 = dx(alf1, alf3, lamb, m, x + h*k3, round(y + h, rnd), round(z + h, rnd))
        pend1 = (k1 + 2*k2 + 2*k3 + k4)/6
        x_temp = x + h*pend1

        k1 = dy(alf0, r0, y, z)
        k2 = dy(alf0, r0, y + h*k1/2, round(z + (h/2), rnd + 1))
        k3 = dy(alf0, r0, y + h*k2/2, round(z + (h/2), rnd + 1))
        k4 = dy(alf0, r0, y + h*k3, round(z + h, rnd))
        pend2 = (k1 + 2*k2 + 2*k3 + k4)/6
        y_temp = y + h*pend2

        k1 = dz(ganm, b, x, y, z)
        k2 = dz(ganm, b, round(x + (h/2), rnd + 1), round(y + (h/2), rnd + 1), z + h*k1/2)
        k3 = dz(ganm, b, round(x + (h/2), rnd + 1), round(y + (h/2), rnd + 1), z + h*k2/2)
        k4 = dz(ganm, b, round(x + h, rnd), round(y + h, rnd), z + h*k3)
        pend3 = (k1 + 2*k2 + 2*k3 + k4)/6
        z_temp = z + h*pend3

        t = round(t + h,rnd)
        x = x_temp
        y = y_temp
        z = z_temp
        ts.append(t)
        xs.append(x)
        ys.append(y)
        zs.append(z)
    return (ts, xs ,ys, zs)

if __name__ == '__main__':
    t = 0
    tn = 150
    x = 1
    y = .1
    z = .001
    h = .01
    alf1 = .1
    alf2 = 1
    alf3 = 0
    b = 1
    ganm = 1
    lamb = .1
    m = .01
    r = .1
    rk = RungeKutta(dx, dy, dz, alf1, alf2, alf3, b, ganm, lamb, m, r, x, y, z, h, t, tn, 3)
    
    ts = rk[0]
    xs = rk[1]
    ys = rk[2]
    zs = rk[3]
    plt.plot(ts, xs, label="Células inactivas")
    plt.plot(ts, ys, label="Células proliferantes")
    plt.plot(ts, zs, label="Anticuerpos")
    plt.legend()
    plt.xlabel("Tiempo")
    plt.ylabel("Densidad Tumoral")
    plt.grid()
    plt.show()