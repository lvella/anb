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
    obst = Matrix(False, lx, ly)
    obst_file = open('anb.obs', 'r')
    for line in obst_file:
        (x, y) = [int(n) - 1 for n in line.split()]
        if 0 <= x < lx and 0 <= y < ly:
            obst[x, y] = True
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

    for y in obst.dims[1]:
        if (not obst[0, y]
            and node[2, 0, y] - t_1 > 0
            and node[5, 0, y] - t_2 > 0
            and node[6, 0, y] - t_2 > 0
            ):

            node[0, 0, y] += t_1
            node[2, 0, y] -= t_1
            node[4, 0, y] += t_2
            node[5, 0, y] -= t_2
            node[6, 0, y] -= t_2
            node[7, 0, y] += t_2

def propagate():
    pass

def bouceback(obst, node, n_hlp):
    pass

def relaxation():
    pass

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

    for time in xrange(1, t_max+1):
        if time % (t_max//10) == 0:
            check_density(node, time)

        

if __name__ == '__main__':
    main()
