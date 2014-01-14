import sys
import operator
import itertools

class Matrix2(list):
    def __init__(self, default, *dims):
        self.dims = dims
        coefs = []

        # Coeficients for calculating linear index
        size = 1
        for d in dims:
            size *= d
            coefs.append(size)
        self.coefs = coefs

        self.extend([default] * size)

    def _key_map(self, key):
        # TODO: to be continued
        return y * self.lx + x

    def __getitem__(self, key):
        return list.__getitem__(self, self._key_map(key))

    def __setitem__(self, key, value):
        return list.__setitem__(self, self._key_map(key), value)

    def __iter__(self):
        iterator = list.__iter__(self)
        return (itertools.islice(iterator, self.lx) for 
                i in xrange(self.ly))

def read_parametrs():
    par_file = open('anb.par', 'r')
    ret = []
    for line in par_file:
        try:
            ret.append(float(line))
        except ValueError:
            pass
    return ret

def read_obstacles(lx, ly):
    obst = Matrix2(False, lx, ly)
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

    #node = # TODO: to be continued

def main():
    if len(sys.argv) >= 3:
        lx = int(sys.argv[1])
        ly = int(sys.argv[2])
    else:
        lx = 30
        ly = 20

    t_max, density, accel, omega, r_rey = read_parametrs()
    obst = read_obstacles(lx, ly)
    for row in obst:
        for cell in row:
            sys.stdout.write('X' if cell else ' ')
        sys.stdout.write('\n')

if __name__ == '__main__':
    main()
