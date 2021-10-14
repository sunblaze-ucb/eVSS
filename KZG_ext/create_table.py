import os

s_time = []
p_time = []
v_time = []
proof_size = []
N = 20

for i in range(5, N):
    
    os.system('./KZG_std ' + str(i) + ' ' + str(i + 1) + ' >> log' + str(i) + '.txt')
    f = open('log' + str(i) + '.txt')

    lines = f.readlines()
    print("Dealing for N = " + str(2**i))
    print(lines)
    p = str(float(lines[-3].split(' ')[-1]))
    v = str(float(lines[-1].split(' ')[-1]))
    s = str(float(lines[1].split(' ')[-1]))
    print(p, v, s)

    p_time.append(p)
    v_time.append(v)
    s_time.append(s)

    print('\n\n\n')


lines = ['', '', '']


for i in range(10, N):
    lines[0] = lines[0] + ' ' + str(p_time[i - 10])
    lines[1] = lines[1] + ' ' + str(v_time[i - 10])
    lines[2] = lines[2] + ' ' + str(s_time[i - 10])

print(lines[0])
print(lines[1])
print(lines[2])

