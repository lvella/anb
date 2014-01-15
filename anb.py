import sys
import operator
import itertools
import functools

def slicer(iterator, size, count):
     return (itertools.islice(iterator, size) for i in xrange(count))

class Matrix(list):
    def __init__(self, default, *dims):
        self.dims = dims

        # Coeficients for calculating linear index
        coefs = []

        size = 1
        for d in dims:
            coefs.append(size)
            size *= d
        self.coefs = coefs
        self.extend([default] * size)

        # Used on iter computation
        sizes = []
        partial = 1
        for d in itertools.islice(reversed(dims), len(dims)-1):
            partial *= d
            sizes.append(partial)
        sizes.reverse()
        self.sizes = sizes
        print sizes

    def _key_map(self, key):
        return sum((a*b for (a,b) in itertools.izip(self.coefs, key)))

    def __getitem__(self, key):
        return list.__getitem__(self, self._key_map(key))

    def __setitem__(self, key, value):
        return list.__setitem__(self, self._key_map(key), value)

    def __iter__(self):
        iterator = list.__iter__(self)
        for i in xrange(len(self.dims) - 1):
            iterator = slicer(iterator, self.dims[i], self.sizes[i])
        return iterator

    def as_list(self):
        return super(Matrix, self)

def read_parametrs():
    par_file = open('anb.par', 'r')
    ret = [f(par_file.readline()) for f in [int, float, float, float, int]]
    return ret

def read_obstacles(lx, ly):
    obst = set()
    obst_file = open('anb.obs', 'r')
    for line in obst_file:
        (x, y) = [int(n) - 1 for n in line.split()]
        if 0 <= x < lx and 0 <= y < ly:
            obst.add((x, y))
        else:
            print '!!! Obstacle out of range, skipped'
            print '!!! lx =', x, ' , ly =', y
            print '!!!'
    return obst

def init_density(lx, ly, density):
    t_0 = density * 4.0 / 9.0
    t_1 = density / 9.0
    t_2 = density / 36.0

    node = Matrix(0.0, 9, lx, ly)
    node[:] = ([t_0] + [t_1] * 4 + [t_2] * 4) * (lx * ly)

    return node

def check_density(node, time):
    n_sum = sum(node.as_list().__iter__())
    print '*** Iteration number = ', time
    print '*** Integral density = ', n_sum
    print '***'

def redistribute(obst, node, accel, density):
    t_1 = density * accel / 9.0
    t_2 = density * accel / 36.0

    for y in xrange(node.dims[2]):
        if ((0, y) not in obst
            and node[3, 0, y] - t_1 > 0
            and node[6, 0, y] - t_2 > 0
            and node[7, 0, y] - t_2 > 0
            ):

            node[1, 0, y] += t_1
            node[3, 0, y] -= t_1
            node[5, 0, y] += t_2
            node[6, 0, y] -= t_2
            node[7, 0, y] -= t_2
            node[8, 0, y] += t_2

def propagate(node, n_hlp):
    lx = node.dims[1]
    ly = node.dims[2]

    for x in xrange(lx):
        x_e = (x+1) % lx
        x_w = (lx + x - 1) % lx

        for y in xrange(ly):
            y_n = (y+1) % ly
            y_s = (ly + x - 1) % ly

            for i, dx, dy in (
                (0, x, y),
                (1, x_e, y),
                (2, x, y_n),
                (3, x_w, y),
                (4, x, y_s),
                (5, x_e, y_n),
                (6, x_w, y_n),
                (7, x_w, y_s),
                (8, x_e, y_s)
            ):
                n_hlp[i, dx, dy] = node[i, x, y]

def bounceback(obst, node, n_hlp):
    for (x, y) in obst:
        for (a, b) in bounceback.permute:
            node[a, x, y] = n_hlp[b, x, y]
bounceback.permute = ((1,3), (2,4), (3,1), (4,2), (5,7), (6,8), (7,5), (8,6))

def relaxation(density, omega, node, n_hlp, obst):
    t_0 = 4.0 / 9.0
    t_1 = 1.0 / 9.0
    t_2 = 1.0 / 36.0
    c_squ = 1.0 / 3.0

    for x in xrange(node.dims[1]):
        for y in xrange(node.dims[2]):
            if (x, y) in obst:
                continue

            d_loc = 0.0
            for i in xrange(9):
                d_loc += n_hlp[i, x, y]

            u_x = (n_hlp[1,x,y] + n_hlp[5,x,y] + n_hlp[8,x,y]
                    - (n_hlp[3,x,y] + n_hlp[6,x,y] + n_hlp[7,x,y])) / d_loc

            u_y = (n_hlp[2,x,y] + n_hlp[5,x,y] + n_hlp[6,x,y]
                    -(n_hlp[4,x,y] + n_hlp[7,x,y] + n_hlp[8,x,y])) / d_loc

            u_squ = u_x * u_x + u_y * u_y

            u_n = [0.0] * 9
            u_n[1] =   u_x
            u_n[2] =         u_y
            u_n[3] = - u_x
            u_n[4] =       - u_y
            u_n[5] =   u_x + u_y
            u_n[6] = - u_x + u_y
            u_n[7] = - u_x - u_y
            u_n[8] =   u_x - u_y

            n_equ = [0.0] * 9
            n_equ[0] = t_0 * d_loc * (1.0 - u_squ / (2.0 * c_squ))
            for i, t in zip(xrange(1,9), [t_1] * 4 + [t_2] * 4):
                n_equ[i] = t * d_loc * (1.0 + u_n[i] / c_squ 
                        + (u_n[i] * u_n[i]) / (2.0 * c_squ * c_squ) 
                        - u_squ / (2.0 * c_squ))

            for i in xrange(9):
                node[i,x,y] = n_hlp[i,x,y] + omega * (n_equ[i] - n_hlp[i,x,y])

def write_velocity():
    pass

def write_results():
    pass

def comp_rey():
    pass

def main():
    if len(sys.argv) >= 3:
        lx = int(sys.argv[1])
        ly = int(sys.argv[2])
    else:
        lx = 30
        ly = 20

    t_max, density, accel, omega, r_rey = read_parametrs()
    obst = read_obstacles(lx, ly)
#    for row in obst:
#        for cell in row:
#            sys.stdout.write('X' if cell else ' ')
#        sys.stdout.write('\n')

    node = init_density(lx, ly, density)
#    for level in node:
#        for row in level:
#            for cell in row:
#                sys.stdout.write('{} '.format(cell))
#            sys.stdout.write('\n')
#        sys.stdout.write('\n===================================\n\n')

    n_hlp = Matrix(0.0, 9, lx, ly)

    # Main loop:
    for time in xrange(1, t_max+1):
        if time % (t_max//10) == 0:
            check_density(node, time)


if __name__ == '__main__':
    main()
