"""
RKF45 HPC version to be run on the Legion cluster
--------------------------------------------------
Changes:
- The state and the Fisher information will be computed at the same time.
"""
from numpy import *
from qutip import *
import matplotlib.pyplot as plt
import time


def RKF45(f, t0, tf, initial_statevec, tol, hmax, hmin, N):
    # f = self.Lindblad
    # a = start time = 0
    # b = end time = self.time
    # x0 = statevec with initial states
    # tol = self.tolerance
    # hmax = self.hmax
    # hmin = self.hmin
    

    t = t0
    x = initial_statevec
    h = hmax
    tf = tf - t0

    # Initialise arrays to be returned - put in initial value
    T = array([t])
    fock_fisher = []
    homodyne_fisher = []
    heterodyne_fisher = []

    # Generate the eigenstates for the various Fisher informations
    fock_eigenstates = []
    homodyne_eigenstates = []
    heterodyne_eigenstates = []

    # Initialise arrays
    fock_probs = []
    fock_diffs = []
    homodyne_probs = []
    homodyne_diffs = []
    heterodyne_probs = []
    heterodyne_diffs = []

    # Define sum variables
    fock_sum = 0
    homodyne_sum = 0
    heterodyne_sum = 0

    # Testing Quadratures
    x_quad = array([])
    p_quad = array([])

    # Fock eigenstates
    #for n in range(0, N):
    #    fock_eigenstates.append((tensor(qeye(N), ket2dm(fock(N, n)))).full())

    # Position eigenstates
    position = (destroy(N) + destroy(N).dag())/sqrt(2.)
    position_eigenstates = position.eigenstates()[1]

    for eigenstate in position_eigenstates:
        homodyne_eigenstates.append((tensor(ket2dm(eigenstate), qeye(N))).full())



    # Momentum eigenstates
    momentum = -1.j*(destroy(N) - destroy(N).dag())/sqrt(2.)
    momentum_eigenstates = momentum.eigenstates()[1]
    for eigenstate in momentum_eigenstates:
        heterodyne_eigenstates.append((tensor(ket2dm(eigenstate), qeye(N))).full())

    #position_quad = (tensor(position, qeye(N))).full()
    #momentum_quad = (tensor(momentum, qeye(N))).full()
    position_quad = (tensor(qeye(N), position)).full()
    momentum_quad = (tensor(qeye(N), momentum)).full()


    """
    Perform an initial calculation of the quads and the Fisher information
    """
    x_quad = append(x_quad, real(trace(dot(x[0],position_quad))))
    p_quad = append(p_quad, real(trace(dot(x[0],momentum_quad))))

    for i in range(0,N):
        #fock_probs.append(real(trace(dot(fock_eigenstates[i],x[0]))))
        #fock_diffs.append(real(trace(dot(fock_eigenstates[i],x[1]))))
        homodyne_probs.append(real(trace(dot(homodyne_eigenstates[i],x[0]))))       
        homodyne_diffs.append(real(trace(dot(homodyne_eigenstates[i],x[1]))))
        heterodyne_probs.append(real(trace(dot(heterodyne_eigenstates[i],x[0]))))
        heterodyne_diffs.append(real(trace(dot(heterodyne_eigenstates[i],x[1]))))


    for i in range(0, N):
        #fock_sum = fock_sum + fock_diffs[i]**2/fock_probs[i]
        homodyne_sum = homodyne_sum + homodyne_diffs[i]**2/homodyne_probs[i]
        heterodyne_sum = heterodyne_sum + heterodyne_diffs[i]**2/heterodyne_probs[i]

    #fock_fisher.append(fock_sum)
    homodyne_fisher.append(homodyne_sum)
    heterodyne_fisher.append(heterodyne_sum)
    #fock_probs = []
    #fock_diffs = []
    homodyne_probs = []
    homodyne_diffs = []
    heterodyne_probs = []
    heterodyne_diffs = []
    #fock_sum = 0
    homodyne_sum = 0
    heterodyne_sum = 0



    # Coefficients in the RK method
    b21 =   2.500000000000000e-01  #  1/4
    b31 =   9.375000000000000e-02  #  3/32
    b32 =   2.812500000000000e-01  #  9/32
    b41 =   8.793809740555303e-01  #  1932/2197
    b42 =  -3.277196176604461e+00  # -7200/2197
    b43 =   3.320892125625853e+00  #  7296/2197
    b51 =   2.032407407407407e+00  #  439/216
    b52 =  -8.000000000000000e+00  # -8
    b53 =   7.173489278752436e+00  #  3680/513
    b54 =  -2.058966861598441e-01  # -845/4104
    b61 =  -2.962962962962963e-01  # -8/27
    b62 =   2.000000000000000e+00  #  2
    b63 =  -1.381676413255361e+00  # -3544/2565
    b64 =   4.529727095516569e-01  #  1859/4104
    b65 =  -2.750000000000000e-01  # -11/40

    # Coefficients used to compute local truncation error estimate.  These
    # come from subtracting a 4th order RK estimate from a 5th order RK
    # estimate.

    r1  =   2.777777777777778e-03  #  1/360
    r3  =  -2.994152046783626e-02  # -128/4275
    r4  =  -2.919989367357789e-02  # -2197/75240
    r5  =   2.000000000000000e-02  #  1/50
    r6  =   3.636363636363636e-02  #  2/55

    # Coefficients used to compute 4th order RK estimate
    c1  =   1.157407407407407e-01  #  25/216
    c3  =   5.489278752436647e-01  #  1408/2565
    c4  =   5.353313840155945e-01  #  2197/4104
    c5  =  -2.000000000000000e-01  # -1/5


    while t < tf:
        if t + h > tf: # Adapt the step size for the final step to make sure we do not overshoot
            h = tf - t

        print t

        # Compute values to compute the truncation error estimate and the 4th order RK estimate.

        k1 = h * f( x )
        k2 = h * f( x + b21 * k1 )
        k3 = h * f( x + b31 * k1 + b32 * k2 )
        k4 = h * f( x + b41 * k1 + b42 * k2 + b43 * k3 )
        k5 = h * f( x + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4 )
        k6 = h * f( x + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5 )

        # Compute the estimate of the local truncation error.  If it's small
        # enough then we accept this step and save the 4th order estimate.
        r = abs( r1 * k1 + r3 * k3 + r4 * k4 + r5 * k5 + r6 * k6 ) / h #Calculate the difference between RK5 and RK4
        if len( shape( r ) ) > 0:
            r = amax( r ) # Find the largest value, this is the largest error
        if r <= tol: # If the largest value is smaller than our tolerance, accept the values
            t = t + h
            x = x + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
            T = append( T, t )
            """
            Test Quadratures
            """
            x_quad = append(x_quad, real(trace(dot(x[0],position_quad))))
            p_quad = append(p_quad, real(trace(dot(x[0],momentum_quad))))

            """
            This is where we calculate the Fisher Information.
            """
            # Fock measurement
            for i in range(0,N):
                #fock_probs.append(real(trace(dot(fock_eigenstates[i],x[0]))))
                #fock_diffs.append(real(trace(dot(fock_eigenstates[i],x[1]))))
                homodyne_probs.append(real(trace(dot(homodyne_eigenstates[i],x[0]))))
                homodyne_diffs.append(real(trace(dot(homodyne_eigenstates[i],x[1]))))
                heterodyne_probs.append(real(trace(dot(heterodyne_eigenstates[i],x[0]))))
                heterodyne_diffs.append(real(trace(dot(heterodyne_eigenstates[i],x[1]))))


            #print fock_diffs
            for i in range(0, N):
                #fock_sum = fock_sum + fock_diffs[i]**2/fock_probs[i]
                homodyne_sum = homodyne_sum + homodyne_diffs[i]**2/homodyne_probs[i]
                heterodyne_sum = heterodyne_sum + heterodyne_diffs[i]**2/heterodyne_probs[i]

            #fock_fisher.append(fock_sum)
            homodyne_fisher.append(homodyne_sum)
            heterodyne_fisher.append(heterodyne_sum)
            #fock_probs = []
            #fock_diffs = []
            homodyne_probs = []
            homodyne_diffs = []
            heterodyne_probs = []
            heterodyne_diffs = []
            #fock_sum = 0
            homodyne_sum = 0
            heterodyne_sum = 0

            """
            Here we calculate the entropy and such if needed.
            """

        # Now compute next step size, and make sure that it is not too big or
        # too small.
        h = h * min( max( 0.84 * ( tol / r )**0.25, 0.1 ), 4.0 )

        if h > hmax:
            h = hmax
        elif h < hmin:
            print "Error: stepsize should be smaller than %e." % hmin
            break
    #End of loop


    # Return the full array at the end of the loop
    return [T, asarray(fock_fisher), asarray(homodyne_fisher), asarray(heterodyne_fisher), x_quad, p_quad]
