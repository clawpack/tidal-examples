import numpy as np
import matplotlib.pyplot as plt
time = np.genfromtxt ('hours', delimiter=",")

randall1102 = np.genfromtxt ('randall1102', delimiter=",")
randall1187 = np.genfromtxt ('randall1187', delimiter=",")

simple1102 =np.genfromtxt ('simple1102', delimiter=",")
simple1187 = np.genfromtxt ('simple1187', delimiter=",")

riemann1102 =np.genfromtxt ('riemann1102', delimiter=",")
riemann1187 =np.genfromtxt ('riemann1187', delimiter=",")

plt.figure(figsize=(13,5))
plt.plot(time,randall1102, label='Randall1102')
plt.plot(time,simple1102, label='simple1102')
plt.plot(time,riemann1102, label='Riemann1102')


plt.xlim(0,72)
plt.ylim(-0.5,0.5)
plt.legend(loc='upper left', fontsize=8)
plt.xlabel('hours')
plt.ylabel('m')
plt.title('NOAA minus GeoClaw at 1102, Dec. 21-23, 2015')
plt.grid(True)
plt.savefig('1102error.png')

plt.figure(figsize=(13,5))

plt.plot(time,randall1187, label='Randall1187')
plt.plot(time,simple1187, label='simple1187')
plt.plot(time,riemann1187, label='Riemann1187')



plt.xlim(0,72)
plt.ylim(-0.5,0.5)
plt.legend(loc='upper left', fontsize=8)
plt.xlabel('hours')
plt.ylabel('m')
plt.title('NOAA minus GeoClaw at 1187, Dec. 21-23, 2015')
plt.grid(True)

plt.savefig('1187error.png')

print("Randall, 1 norm error at 1102:", np.linalg.norm(randall1102,1))
print("Simple, 1 norm error at 1102:",np.linalg.norm(simple1102,1))
print("Riemann, 1 norm error at 1102:",np.linalg.norm(riemann1102,1))

print("Randall, 2 norm error at 1102:", np.linalg.norm(randall1102,2))
print("Simple, 2 norm error at 1102:",np.linalg.norm(simple1102,2))
print("Riemann, 2 norm error at 1102:",np.linalg.norm(riemann1102,2))
