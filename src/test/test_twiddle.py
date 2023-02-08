import math
import cmath
from sys import argv

def twiddle(forward, k, N):
    return cmath.exp(2j * k * (math.pi / N)) if forward else cmath.exp(-2j * k * (math.pi / N))


if __name__ == "__main__":
    N, forward = map(lambda x: int(x), argv[1:])
    print(N, forward)
    for i in range(int(N / 2)):
        twid = twiddle(forward, i, N)
        print("{:.5f} {:.5f}".format(twid.real, twid.imag))
