mod = 21888242871839275222246405745257275088548364400416034343698204186575808495617
import random
def fastpow(x, y, mod):
    ret = 1
    while(y != 0):
        if(y % 2 == 1):
            ret = (ret * x) % mod
        x = (x * x) % mod
        y = y // 2
    return ret
def find_root_of_unity(mod, n):
    while(True):
        x = random.randint(0, mod - 1)
        g = fastpow(x, (mod - 1) // n, mod)
        if(fastpow(g, n // 2, mod) != 1):
            print(g)
            return g
        print('fail')
g = find_root_of_unity(mod, 2**28)
if(fastpow(g, 2 ** 27, mod) != 1 and fastpow(g, 2 ** 28, mod) == 1):
    print('test succ')
    for i in range(29):
        print(fastpow(g, 2 ** i, mod), i)
else:
    print('test fail')