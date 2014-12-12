from myhdl import * 
import numpy as np 
from math import e, pi, log
import matplotlib.pyplot as plt

t_state = enum('INIT', 'DATA_IN', 'COMPUTE', 'COMPUTE_INDEX', 'COMPUTE_MULT', 'DATA_OUT')

#########################CHANGES NEEDED IF N!=8###############################
def FFT(clk, reset, start, data_valid,  
        #ADD PORTS HERE
        #In_real
        xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7,
        #In_imag
        xi0, xi1, xi2, xi3, xi4, xi5, xi6, xi7,
        #Out_real
        zr0, zr1, zr2, zr3, zr4, zr5, zr6, zr7,
        #Out_imag
        zi0, zi1, zi2, zi3, zi4, zi5, zi6, zi7, 
        N_points, Q):
##############################################################################


    points = N_points
    half_points = points>>1
    logpoints = int(log(points,2))
    n_bits = len(xr0)

    #Declaration of the buffers used to go through the Butterfly
    xr_buf = [Signal(modbv(0, min=-2**(n_bits-1), max=2**(n_bits-1))) for i in
            range(points)]
    xi_buf = [Signal(modbv(0, min=-2**(n_bits-1), max=2**(n_bits-1))) for i in
            range(points)]


    state = Signal(t_state.INIT)
    #Level of FFT --> total number of levels = log(points)
    level = Signal(intbv(0, min=0, max= logpoints+1))
    #Step to get the right pair, depends on the level
    step = Signal(intbv(1, min=1, max=2**logpoints+1))
    compute_index = Signal(intbv(0, min=0, max=points))
    fft_index = Signal(intbv(0, min=0, max= half_points+1))
    xr_buf_idx = Signal(intbv(0)[5:])
    increm = Signal(intbv(0)[7:])
    #Signals for twiddles
    ur = Signal(modbv(0,min=-2**(31), max=2**(31)))
    ui = Signal(modbv(0,min=-2**(31), max=2**(31)))
    #Signals for products
    prod01 = Signal(modbv(0,min=-2**(31), max=2**(31)))
    prod02 = Signal(modbv(0,min=-2**(31), max=2**(31)))
    prod11 = Signal(modbv(0,min=-2**(31), max=2**(31)))
    prod12 = Signal(modbv(0,min=-2**(31), max=2**(31)))
    #Used to get Twiddles in the W arrays transformed in tuples
    level_tuple = Signal(intbv(0, min=0, max= (logpoints)*points+1))
    

#####
    #Prepare the twiddles as follows (Example for N=8): 
    #        0     0     0     0     
    # W = [ W     W     W     W   
    #        8     8     8     8
    #        0     0     2     2     
    #       W     W     W     W   
    #        8     8     8     8
    #        0     1     2     3     
    #       W     W     W     W   ] 
    #        8     8     8     8

#########################CHANGES NEEDED IF N!=8###############################
    #For 8 points: 
    tw_index = [[0,0,0,0],
                [0,0,2,2],
                [0,1,2,3]]    
##############################################################################

    W = np.ones([logpoints,points>>1])+np.ones([logpoints,points>>1])*1j
    Wr =  np.zeros([logpoints, points>>1], np.int32)
    Wi =  np.zeros([logpoints, points>>1], np.int32)
    Tw =  np.zeros([logpoints, points>>1], np.int32)
    for k in range(logpoints):
        #index=modbv(0, min=0, max=points>>1)
        for i in range(points>>1): 
            index = tw_index[k][i]
            W[k][i] = e**(-2j*pi*int(index)/points)
            #Tw[k][i] = index
            #index+=points>>(k+1)
            Wr[k][i] = W[k][i].real * 2**Q 
            Wi[k][i] = W[k][i].imag * 2**Q
    Wr = tuple(map(int, Wr.reshape((1,logpoints*(points>>1)))[0]))
    Wi = tuple(map(int, Wi.reshape((1,logpoints*(points>>1)))[0]))
##### 


    @always_seq(clk.posedge, reset)
    def compute(): 
        if state == t_state.INIT:
            data_valid.next = 0
            if start == 1: 
                state.next = t_state.DATA_IN
      
        #fill the buffers in the correct input order of the Butterfly: 
        #Example (N=8): xr[1] --> xr[4]  (001 --> 100)
        elif state == t_state.DATA_IN: 
            state.next = t_state.COMPUTE_INDEX
#########################CHANGES NEEDED IF N!=8###############################
            xr_buf[0].next = xr0
            xr_buf[4].next = xr1
            xr_buf[2].next = xr2
            xr_buf[6].next = xr3
            xr_buf[1].next = xr4
            xr_buf[5].next = xr5
            xr_buf[3].next = xr6
            xr_buf[7].next = xr7
            xi_buf[0].next = xi0
            xi_buf[4].next = xi1
            xi_buf[2].next = xi2
            xi_buf[6].next = xi3
            xi_buf[1].next = xi4
            xi_buf[5].next = xi5
            xi_buf[3].next = xi6
            xi_buf[7].next = xi7
##############################################################################
            

        #To avoid a critical path exceding the timing constraints, 3 states
        #are used to execute the FFT 

        #State1 : Prepare the indeces
        elif state == t_state.COMPUTE_INDEX:
            increm.next = step+step
            if level < level.max-1: 
                #print('level %d of %d'%(level, level.max-2))
                #print('step: %d' %step)
                #print('increm: %d' %increm)
                if fft_index < fft_index.max-1: 
                    ur.next = Wr[level_tuple+fft_index] 
                    ui.next = Wi[level_tuple+fft_index]
                    xr_buf_idx.next =  compute_index+step
                    state.next =t_state.COMPUTE_MULT
                else:  
                    compute_index.next = 0
                    level.next=level+1
                    fft_index.next=0
                    level_tuple.next=level_tuple+half_points
                    step.next=step+step
                    #print('------------NEXT LEVEL--------------') 
            else: 
                state.next = t_state.DATA_OUT
 
        #State2 : Compute the products   
        elif state == t_state.COMPUTE_MULT:
            prod01.next = xr_buf[xr_buf_idx]*ur
            prod02.next = xi_buf[xr_buf_idx]*ui
            prod11.next = xr_buf[xr_buf_idx]*ui
            prod12.next = xi_buf[xr_buf_idx]*ur
            state.next = t_state.COMPUTE
        
        #State3 : Compute the new FFT value
        elif state == t_state.COMPUTE:
            #print('W = %d + i(%d)'%(ur, ui))
            #print('computing: x[%d] & x[%d]' %(compute_index, compute_index+step)) 
            
            prod0 = modbv(prod01 - prod02, min=-2**31, max=2**31)
            prod1 = modbv(prod11 + prod12, min=-2**31, max=2**31)
            xr_buf[compute_index].next = xr_buf[compute_index] + prod0[32:16]
            xi_buf[compute_index].next = xi_buf[compute_index] + prod1[32:16]
            xr_buf[compute_index+step].next = xr_buf[compute_index] - \
                prod0[32:16]
            xi_buf[compute_index+step].next = xi_buf[compute_index] - \
                prod1[32:16]

            #print('xr[%d] = %d'%(compute_index, xr_buf[compute_index]+prod0))
            #print('xi[%d] = %d'%(compute_index, xi_buf[compute_index]+prod1))
            #print('xr[%d] = %d'%(compute_index+step,
            #    xr_buf[compute_index+step]-prod0))
            #print('xi[%d] = %d'%(compute_index+step,
            #    xi_buf[compute_index+step]-prod1))

            compute_index.next = (compute_index+increm)%(points-1)
            fft_index.next=fft_index+1
            state.next = t_state.COMPUTE_INDEX


        #Assign the buffers to the outputs
        elif state == t_state.DATA_OUT:
            data_valid.next = 1 
#########################CHANGES NEEDED IF N!=8###############################
            zr0.next = xr_buf[0]
            zr1.next = xr_buf[1]
            zr2.next = xr_buf[2]
            zr3.next = xr_buf[3]
            zr4.next = xr_buf[4]
            zr5.next = xr_buf[5]
            zr6.next = xr_buf[6]
            zr7.next = xr_buf[7]
            zi0.next = xi_buf[0]
            zi1.next = xi_buf[1]
            zi2.next = xi_buf[2]
            zi3.next = xi_buf[3]
            zi4.next = xi_buf[4]
            zi5.next = xi_buf[5]
            zi6.next = xi_buf[6]
            zi7.next = xi_buf[7]
##############################################################################
            level.next = 0
            level_tuple.next = 0
            step.next = 1
            state.next = t_state.INIT

    return compute 

def compile_FFT(): 
    n_bits=16
    Q=16
    clk = Signal(bool(0))
    reset = ResetSignal(0, active=1, async=True)
    start = Signal(bool(0))
    data_valid = Signal(bool(0))
#########################CHANGES NEEDED IF N!=8###############################
    N_points=8
    xr0 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr1 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr2 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr3 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr4 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr5 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr6 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xr7 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi0 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi1 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi2 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi3 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi4 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi5 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi6 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    xi7 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr0 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr1 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr2 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr3 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr4 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr5 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr6 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zr7 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi0 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi1 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi2 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi3 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi4 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi5 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi6 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
    zi7 = Signal(intbv(0, min=-2**(n_bits-1), max=2**(n_bits-1)))
       
    toVHDL(FFT, clk, reset, start, data_valid, 
            xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7,
            xi0, xi1, xi2, xi3, xi4, xi5, xi6, xi7,
            zr0, zr1, zr2, zr3, zr4, zr5, zr6, zr7,
            zi0, zi1, zi2, zi3, zi4, zi5, zi6, zi7,
            N_points, Q) 
##############################################################################



#SIMULATION
def FFT_tb():
    
    HALF_PERIOD = delay(5)
    n_bits=16
    Q=16
    
    clk = Signal(bool(0))
    reset = ResetSignal(0, active=1, async=True)
    start = Signal(bool(0))
    data_valid = Signal(bool(0))
    
#########################CHANGES NEEDED IF N!=8###############################
    N_points=8
   
    [xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7,
    xi0, xi1, xi2, xi3, xi4, xi5, xi6, xi7,
    zr0, zr1, zr2, zr3, zr4, zr5, zr6, zr7,
    zi0, zi1, zi2, zi3, zi4, zi5, zi6, zi7] = [Signal(intbv(0,
        min=-2**(n_bits-1), max=2**(n_bits-1))) for i in range(N_points*4)]
    
    #Can ONLY be usedto simplify the simulation
    fft_in_bus = [xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7]
    fft_re_bus = [zr0, zr1, zr2, zr3, zr4, zr5, zr6, zr7]
    fft_im_bus = [zi0, zi1, zi2, zi3, zi4, zi5, zi6, zi7]

    DUT = FFT(clk, reset, start, data_valid,
                xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7,
                xi0, xi1, xi2, xi3, xi4, xi5, xi6, xi7,
                zr0, zr1, zr2, zr3, zr4, zr5, zr6, zr7,
                zi0, zi1, zi2, zi3, zi4, zi5, zi6, zi7,
                N_points, Q)
##############################################################################

    @always(HALF_PERIOD)
    def clockGen():
        clk.next = not clk

    raw = [-3, -2, -1,  0,  0,  1,  2,  3]   


    @instance
    def tb_stim():
        reset.next = True
        yield clk.posedge
        reset.next = False
        i = 0
        for sig in fft_in_bus: 
            sig.next = raw[i] 
            i+=1
        yield clk.posedge
        yield clk.negedge
        yield clk.negedge
        start.next = 1
        yield clk.negedge
        start.next = 0
        for i in range(50):
            yield clk.negedge
        
        X = np.zeros(len(fft_im_bus), dtype='complex')
        x = np.zeros(len(fft_im_bus))
        for i in range(len(x)):
            X[i] = fft_re_bus[i] + fft_im_bus[i]*1j
            x[i] = fft_in_bus[i]
        
        X_np = np.fft.fftn(x)
        print(X_np)
        plt.plot(np.abs((X)))
        plt.plot(np.abs(X_np))
        plt.show()

        raise StopSimulation

    return DUT, clockGen, tb_stim

if __name__ == "__main__":
    compile_FFT()
    sim = Simulation(traceSignals(FFT_tb))
    sim.run()
