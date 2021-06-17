
def f(x):
    return 0.00021 * x * x - 0.00019 * x - 0.00053

delta = 0
for i in range(11):
    print(i + 10, float(f(2 ** (i + 10)))  , '\\\\')