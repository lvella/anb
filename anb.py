import sys

def read_parametrs():
    par_file = open('anb.par', 'r')
    return [float(line) for line in par_file]

def read_obstacles():
    pass # TODO: to be continued

def main():
    if len(sys.argv) >= 3:
        lx = int(sys.argv[1])
        ly = int(sys.argv[2])
    else:
        lx = 30
        ly = 20

    t_max, density, accel, omega, r_rey = read_parametrs()
    obst = read_obstacles(lx, ly)
