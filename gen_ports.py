from math import log
import numpy as np
import matplotlib.pyplot as plt

######
#Change N to the desired number of points
N = 8
######


points = N 
logpoints=int(log(points, 2))


#Input values have to be reorder, i.e. their index inverted: 
#Example (N=8): xr[1] --> xr[4]  (001 --> 100)
n_index = logpoints 
index_rev = [0]*points
for i in range(points): 
    index_rev[i] = sum(1<<(n_index-1-j) for j in range(n_index) if i>>j&1)
index_rev = tuple(index_rev)

print("\nINPUTS and OUTPUTS ports:\n") 
print('def FFT(clk, reset, start, data_valid, ')
print("xr"+", xr".join(map(str, range(N)))+",")
print("xi"+", xi".join(map(str, range(N)))+",")
print("zr"+", zr".join(map(str, range(N)))+",")
print("zi"+", zi".join(map(str, range(N)))+",")
print('N_points, Q):')

print("\nAssign input ports to internal buffer:\n")
for i in range(N):
    print("xr_buf[%d].next = xr%d" %(index_rev[i],i))
for i in range(N):
    print("xi_buf[%d].next = xi%d" %(index_rev[i],i))
    
print("\nAssign internal buffer to output ports :\n")
for i in range(N):
    print("zr%d.next = xr_buf[%d]" %(i,i))
for i in range(N):
    print("zi%d.next = xi_buf[%d]" %(i,i))

print("\nSignal declaration:\n")
for i in range(N):
    print("xr%d = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))" %(i))
for i in range(N):
    print("xi%d = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))" %(i))
for i in range(N):
    print("zr%d = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))" %(i))
for i in range(N):
    print("zi%d = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))" %(i))


