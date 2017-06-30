"""
Solver HPC Version
-------------------------
In this version of the code, we have optimised it to run on the Legion cluster.
We have:
- Removed all output
- Reorganised it so that after each iteration, the saved state is used for calculating all possible Fisher informations.
- Changed it so that states are not saved.
"""

import gc
from qutip import *
import yaml
from numpy import *
import matplotlib.pyplot as plt
from scipy import constants as cp
import os
from matplotlib import rc
import datetime
import time
from scipy.special import factorial
from math import isnan as mathisnan
from multiprocessing import Process, freeze_support, Pool
from RKF45 import RKF45 
import copy_reg
import types
import tables


class Solver(object):
    """
    Description: Class that includes the simulation run and more.  
    """
    def __init__(self, config):
        # Extract arguments from the config file
        self.args = {}
        if type(config) == str:
            with open(config) as cfile:
                self.args.update(yaml.load(cfile))
        elif type(config) == dict:
            self.args.update(config)
        else:
            print "Failed to load config arguments"

        # Assign values to constants
        self.function = str(self.args['function'])
        self.state = str(self.args['state'])
        self.squeezing = float(self.args['squeezing'])
        self.displace = float(self.args['displace'])
        self.start_time = float(self.args['start_time'])
        self.N = int(self.args['N'])
        self.M = int(self.args['M'])
        self.alpha = float(self.args['alpha'])
        self.beta = float(self.args['beta'])
        self.chi = float(self.args['chi'])
        if self.args['time'] == 6.24:
            self.time = 2.*pi
        else:
            self.time = float(self.args['time'])
        self.gbar = float(self.args['gbar'])
        self.k = float(self.args['k'])
        self.tolerance = float(self.args['tolerance'])
        self.folder = str(self.args['folder'])
        self.h = float(self.args['h'])
        self.steps = int(self.time/self.h)
        self.hmax = float(self.args['hmax'])
        self.hmin = float(self.args['hmin'])
        self.gamma = float(self.args['gamma'])

        print "Hilbert space is: " + str([self.N, self.M])

        # Select the function to run
        if self.function == 'state_evolution':
            print "Starting function: state_evolution"
            self.state_evolution()
        if self.function == 'RK4F':
            print "Starting function: RK4F"
            self.run_RKF45()
        if self.function == 'triple_state_evolution':
            print "Starting function: triple_state_evolution"
            self.triple_state_evolution()


    def initialise_coherent(self, arg, objtype = 'matrix'):
        """
        Parameters:
            arg: string, determines type 'ket' or 'dm' of outputted state
        Ouput:  two-system state at t = self.start_time
        Note: State is initialised by creating an array of fock states with various n. The array is then summed and the state outputted. 
        """
        terms = []
        eta = 1. - exp(-1j*self.start_time)
        for n in range(0, self.N):
            terms.append(exp(-absolute(self.alpha)**2/2.)*self.alpha**n/(sqrt(factorial(n)))*exp(1j*(self.start_time-sin(self.start_time))*(self.k*n - self.gbar)**2)*exp(((self.k*n - self.gbar)*eta*exp(1j*self.start_time)*self.beta - (self.k*n - self.gbar)*(1 - exp(1j*self.start_time))*exp(-1j*self.start_time)*self.beta)/2.)*tensor(fock(self.N,n), coherent(self.N, self.beta*exp(-1j*self.start_time) + (self.k*n- self.gbar)*(1 - exp(-1j*self.start_time)))))
        state = sum(terms)
        if arg =='ket':
            if objstype == 'numpy':
                return state/sqrt(state.overlap(state)).full()
            if objtype == 'qobj':
                return state/sqrt(state.overlap(state))
        if arg == 'dm':
            if objtype == 'numpy':
                return (state*state.dag()/state.overlap(state)).full() # turn into dm and normalise.
            if objtype == 'qobj':
                return (state*state.dag()/state.overlap(state)) # turn into dm and normalise.


   
    def initialise_d_start_time(self, arg):
        """
        Initialise the derivative of the state with respect to gbar at t = start_time
        """
        if arg == 'dm':
            alpha_terms = []
            beta_terms = []
            for n in range(0, (self.N)):
                for m in range(0, (self.N)):
                    alpha_terms.append(
                    -exp(- absolute(self.alpha)**2)
                    * ((conj(self.alpha)**m *self.alpha**n)/(sqrt(factorial(n)*factorial(m)))
                    * exp((self.start_time - sin(self.start_time))*1j*(self.k**2*(n**2 - m**2) - 2.*self.k*self.gbar*(n-m)))
                    * (2.*pi*(self.start_time - sin(self.start_time))*1j*self.k*(n-m))*fock(self.N, n)*fock(self.N, m).dag()))
                    beta_terms.append(exp(-absolute(self.beta)**2)*(self.beta**n*conj(self.beta)**m/(sqrt(factorial(n)*factorial(m)))*fock(self.N,n)*fock(self.M, m).dag()))
            drho = tensor(sum(alpha_terms), sum(beta_terms))
            return drho.full()
        if arg == 'ket':
            alpha_terms = []
            beta_terms = []
            for n in range(0, (self.N)):
                alpha_terms.append(exp(-absolute(self.alpha)**2/2.)*conj(self.alpha)**n/(sqrt(factorial(n)))*exp((self.start_time - sin(self.start_time))*1j*(self.k**2 * n**2 - 2.*self.k*self.gbar*n))*((self.start_time - sin(self.start_time))*2.*1j*self.k*n)*fock(self.N,n))
                beta_terms.append(exp(-absolute(self.beta)**2/2.)*self.beta**n/(sqrt(factorial(n)))*fock(self.N,n))
            return tensor(sum(alpha_terms), sum(beta_terms)).full()

    def initialise_squeezed(self, s, arg, objtype = 'matrix'):
        if arg == 'ket':
            if objtype == 'qobj':
                d = displace(self.N, self.displace)
                sq = squeeze(self.N, self.squeezing)
                print sq
                state = tensor(d*s*basis(self.N, 0), d*s*basis(self.N, 0))
                state = state/sqrt(state.overlap(state))
                print state.overlap(state)
                return state

            if objtype == 'matrix':
                return tensor(squeeze(self.N, self.squeezing), squeeze(self.M, self.squeezing)).full()

    def initialise_squeezed2(self, z, arg, objtype = 'qobj'):
        terms = []
        for n in range(0, self.N/2):
            terms.append(((z-1)/(2*(z+1)))**n *sqrt(factorial(2*n))/(factorial(n))*basis(self.N, n))
        state = sum(terms)
        state = state/(sqrt(state.overlap(state)))
        print state.overlap(state)
        return state

    def run_RKF45(self):
        # Create the vector
        states = []
        dstates = []
        statevec = array([self.initialise_coherent('dm'), self.initialise_d_start_time('dm')])
        print "Starting simulation."
        # Set decoherence
        self.Lind = sqrt(self.chi)*tensor(destroy(self.N), qeye(self.M))
        self.LinddagLind = self.Lind.dag()*self.Lind
        self.LinddagLind = self.LinddagLind.full()
        self.Linddag = self.Lind.dag().full()
        self.Lind = self.Lind.full()

        tic = time.clock()

        # Create folder destination, check if exists if not create
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')
        foldername = self.folder + "/simulation" + st
        if not os.path.exists(foldername):
            os.makedirs(foldername)

        # Add the current config file into the simulation folder
        with open(foldername + "/config", 'w') as outfile:
            outfile.write(yaml.dump(self.args, default_flow_style = False))

        results = RKF45(self.Lindblad, self.start_time, self.time, statevec, self.tolerance, self.hmax, self.hmin, self.N)
        print "Simulation completed."
        toc = time.clock()
        print "Time required: " + str((toc-tic)/60) + " minutes"

        # Save results to files
        save(foldername + "/times", results[0])
        save(foldername + "/fock_fisher", results[1])
        save(foldername + "/homodyne_fisher", results[2])
        save(foldername + "/heterodyne_fisher", results[3])
        save(foldername + "/x_quads", results[4])
        save(foldername + "/p_quads", results[5])

        plt.figure(figsize=(10,7.5))
        plt.plot(results[0], results[3])
        plt.show()


    def analytic_U(self, time):
        """
        Parameters: config file and time
        Output: Time evolution operator U generated at time. 
        Note: Used for making sure numerical evolution was giving the correct value. 
        """
        a = tensor(destroy(self.N), qeye(self.N))
        b = tensor(qeye(self.N), destroy(self.N))

        eta = 1 - exp(- 1j*time)
        U = (1j*(self.k*a.dag()*a - self.gbar)**2. *(time - sin(time))).expm()*((self.k*a.dag()*a - self.gbar)*(eta*b.dag() - conj(eta)*b)).expm()*(- 1j*time*b.dag()*b).expm()
        return U


    def state_evolution(self):
        """
        Description: Evolves a state five times and feeds them to a Fisher information function
        Parameters: None
        Output: Saves the resulting Fisher informatio, times and the config to path + self.foldername
        """
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')


        a = tensor(destroy(self.N), qeye(self.N))
        b = tensor(qeye(self.N), destroy(self.N))

        Hamiltonian0 =   b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) 
        times = linspace(self.start_time, self.time, self.steps)
        
        #state0 = tensor(coherent(self.N, self.alpha), coherent(self.N, self.beta))

        if self.state == 'coherent':
            state0 = self.initialise_coherent('ket', 'qobj')
        if self.state == 'squeezed':
            state0 = self.initialise_squeezed2(self.squeezing, 'ket', 'qobj')
        #state0 = tensor(coherent(self.N, 1.), coherent(self.N, 1.))

        g0 = self.gbar

        self.gbar = g0 - 2*self.h
        Hamiltonianm2 = b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) 
        if self.state == 'coherent':
            statem2 = self.initialise_coherent('ket', 'qobj')
        if self.state == 'squeezed':
            statem2 = self.initialise_squeezed2(self.squeezing, 'ket', 'qobj')
        
        self.gbar = g0 - self.h
        Hamiltonianm1 = b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) 
        if self.state == 'coherent':
            statem1 = self.initialise_coherent('ket', 'qobj')
        if self.state == 'squeezed':
            statem1 = self.initialise_squeezed2(self.squeezing, 'ket', 'qobj')


        self.gbar = g0 + self.h
        Hamiltonian1 =  b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) 
        if self.state == 'coherent':
            state1 = self.initialise_coherent('ket', 'qobj')
        if self.state == 'squeezed':
            state1 = self.initialise_squeezed2(self.squeezing, 'ket', 'qobj')

        self.gbar = g0 + 2*self.h
        Hamiltonian2 =  b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) 
        if self.state == 'coherent':
            state2 = self.initialise_coherent('ket', 'qobj')
        if self.state == 'squeezed':
            state2 = self.initialise_squeezed2(self.squeezing, 'ket', 'qobj')


        # Set options
        options = Options()

        # Set decoherence 
        deco = sqrt(self.chi)*tensor(destroy(self.N), qeye(self.N))
        
        psi0 = []
        psi1 = []
        psi2 = []
        psim1 = []
        psim2 = []

        # Perform simulation, trace out the states then delete the results. This is done to preserve memory. 


        print "Starting simulation"
        results0 = mesolve(Hamiltonian0, state0, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        for state in results0.states:
            psi0.append(state.ptrace(0))
        del results0
        gc.collect()

        results1 = mesolve(Hamiltonian1, state1, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        for state in results1.states:
            psi1.append(state.ptrace(0))
        del results1
        gc.collect()

        results2 = mesolve(Hamiltonian2, state2, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        for state in results2.states:
            psi2.append(state.ptrace(0))
        del results2
        gc.collect()


        resultsm1 = mesolve(Hamiltonianm1, statem1, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        for state in resultsm1.states:
            psim1.append(state.ptrace(0))
        del resultsm1
        gc.collect()

        resultsm2 = mesolve(Hamiltonianm2, statem2, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        for state in resultsm2.states:
            psim2.append(state.ptrace(0))
        del resultsm2
        gc.collect()

        # Check quadratures for clarity

        x_exp = []
        p_exp = []
        x = (destroy(self.N) + destroy(self.N).dag())/sqrt(2.)
        p = 1.j * (destroy(self.N).dag() - destroy(self.N))/sqrt(2.)

        for state in psi0:
            x_exp.append(expect(x, state))
            p_exp.append(expect(p, state))
        plt.figure(figsize=(10,7.5))
        plt.plot(x_exp, p_exp)
        plt.show()
 
        tic = time.clock()


        print "Calculating Fisher information"
        [fisher_result_pos, fisher_result_mom] = self.triple_Fisher(psim2, psim1, psi0, psi1, psi2)

        toc = time.clock()
        print "Time required to do that: " + str((toc-tic)/60) + " minutes, or " + str((toc-tic)/(60*60)) + "hours."

        plt.figure(figsize=(10,7.5))
        plt.plot(times, fisher_result_pos)
        plt.show()
        plt.plot(times, fisher_result_mom)
        plt.show()

        # Create folder if it does not exist
        foldername = self.folder + "/simulation" + st
        if not os.path.exists(foldername):
            os.makedirs(foldername)

        with open(foldername + '/config', 'w') as outfile:
            outfile.write(yaml.dump(self.args, default_flow_style=False))
        qsave(psi0, foldername + "/psi0")
        qsave(psi1, foldername + "/psi1")
        qsave(psim1, foldername + "/psim1")
        save(foldername + "/fisher_position", fisher_result_pos)
        save(foldername + "/fisher_momentum", fisher_result_mom)
        save(foldername + "/times", times)


    def triple_state_evolution(self):
        """
        In this section, we let a three-mode quantum state evolve. The last state is a quantum vacuum state that we will measure on. 
        The question is whether we should let the qutip solver handle this one. The problem is that we need the same time steps. 
        So: qutip objects: good solver, can trace out things, can easily take the dagger
        Numpy arrays: Fast solver, ensures same time steps. 
        So the only problem is the times. As long as the solver gives that to us, it's worth keeping the objects. 
        """
        # Take a time-stamp
        st = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d-%H.%M.%S')

        a = tensor(destroy(self.N), qeye(self.N), qeye(self.N))
        b = tensor(qeye(self.N), destroy(self.N), qeye(self.N))
        c = tensor(qeye(self.N), qeye(self.N), destroy(self.N))

        times = linspace(self.start_time, self.time, self.steps)

        # Create Hamiltonian and state for gbar = starting value
        # Add a vacuum state as system C
        Hamiltonian0 =  b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) + self.gamma*(a.dag()*c + a*c.dag())
        state0 = tensor(self.initialise_start_time('ket','qobj'), fock(self.N, 0))

        g0 = self.gbar
        self.gbar = g0 - 2*self.h
        Hamiltonianm2 = b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) + self.gamma*(a.dag()*c + a*c.dag())
        statem2 = tensor(self.initialise_start_time('ket', 'qobj'), fock(self.N,0))

        self.gbar = g0 - self.h
        Hamiltonianm1 = b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) + self.gamma*(a.dag()*c + a*c.dag())
        statem1 = tensor(self.initialise_start_time('ket', 'qobj'), fock(self.N,0))

        self.gbar = g0 + self.h
        Hamiltonian1 =  b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) + self.gamma*(a.dag()*c + a*c.dag())
        state1 = tensor(self.initialise_start_time('ket','qobj'), fock(self.N,0))

        self.gbar = g0 + 2*self.h
        Hamiltonian2 =  b.dag()*b + (b.dag() + b)*(self.gbar - self.k*a.dag()*a) + self.gamma*(a.dag()*c + a*c.dag())
        state2 = tensor(self.initialise_start_time('ket','qobj'), fock(self.N,0))



        psi0 = []
        psi1 = []
        psi2 = []
        psim1 = []
        psim2 = []

        # Add zero decoherence
        deco = []

        # Solve each system for a psi with different g 
        results0 = mesolve(Hamiltonian0, state0, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        print "Tracing out states"
        for state in results0.states:
            psi0.append(state.ptrace(2))
        del results0
        gc.collect()

        results1 = mesolve(Hamiltonian1, state1, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        print "Tracing out states"
        for state in results1.states:
            psi1.append(state.ptrace(2))
        del results1
        gc.collect()

        results2 = mesolve(Hamiltonian2, state2, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        print "Tracing out states"
        for state in results2.states:
            psi2.append(state.ptrace(2))
        del results2
        gc.collect()

        resultsm1 = mesolve(Hamiltonianm1, statem1, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        print "Tracing out states"
        for state in resultsm1.states:
            psim1.append(state.ptrace(2))
        del resultsm1
        gc.collect()


        resultsm2 = mesolve(Hamiltonianm2, statem2, times, c_ops = deco, e_ops = [], args = [], progress_bar = True)
        print "Tracing out states"
        for state in resultsm2.states:
            psim2.append(state.ptrace(2))
        del resultsm2
        gc.collect()

        # Check that we get the correct quadratures
        x = (destroy(self.N) + destroy(self.N).dag())/sqrt(2.)
        p = 1.j * (destroy(self.N).dag() - destroy(self.N))/sqrt(2.)
        x_exp = []
        p_exp = []
        for state in psi0:
            x_exp.append(expect(x, state))
            p_exp.append(expect(p, state))
        plt.figure(figsize=(10,7.5))
        plt.plot(x_exp, p_exp)
        plt.show()
        
        print "Calculating Fisher information by tracing out states"
        tic = time.clock()

        # Here we trace out the states. I wonder however if it wouldn't be faster to just calculate the expectation values. 
        print "Tracing out states"
        # Trace out system AB, to be left with the light state

        fisher_result_pos, fisher_result_mom = self.triple_Fisher(psim2, psim1, psi0, psi1, psi2)

        toc = time.clock()
        print "Time required to do that: " + str((toc-tic)/60) + " minutes, or " + str((toc-tic)/(60*60)) + "hours."

        plt.figure(figsize=(10,7.5))
        plt.plot(times, fisher_result_pos)
        plt.show()
        plt.plot(times, fisher_result_mom)
        plt.show()

        # Create folder if it does not exist
        foldername = self.folder + "/simulation" + st
        if not os.path.exists(foldername):
            os.makedirs(foldername)

        with open(foldername + '/config', 'w') as outfile:
            outfile.write(yaml.dump(self.args, default_flow_style=False))
        qsave(psi0, foldername + "/psi0")
        qsave(psi1, foldername + "/psi1")
        qsave(psim1, foldername + "/psim1")
        save(foldername + "/fisher_leaky_position", fisher_result_pos)
        save(foldername + "/fisher_leaky_momentum", fisher_result_mom)
        save(foldername + "/times", times)



    def triple_Fisher(self, psim2, psim1, psi0, psi1, psi2):
        """
        Description: Calculates the Fisher information for the specified measurement. 
        Parameters:
            psim2: States for a run with self.g - 2*h
            psim1: States for a run with self.g - h
            psi0: States for a run with self.g
            psi1: States for a run with self.g + h
            psi2: States for a run with self.g + 2*h
        Output: An array containig the Fisher information for each timestep
        """
        c_position = (destroy(self.N) + destroy(self.N).dag())/sqrt(2.)
        c_pos_eigen = c_position.eigenstates()
        
        c_momentum = 1j*(destroy(self.N).dag() - destroy(self.N))/sqrt(2.)
        c_mom_eigen = c_momentum.eigenstates()

        probs0_pos = []
        probs0_mom = []
        probs1_pos = []
        probs1_mom = []
        probs2_pos = []
        probs2_mom = []
        probsm1_pos = []
        probsm1_mom = []
        probsm2_pos = []
        probsm2_mom = []


        for state in psim2:
            probsm2_pos.append([expect(ket2dm(x), state) for x in c_pos_eigen[1]])
            probsm2_mom.append([expect(ket2dm(x), state) for x in c_mom_eigen[1]])
        for state in psim1:
            probsm1_pos.append([expect(ket2dm(x), state) for x in c_pos_eigen[1]])
            probsm1_mom.append([expect(ket2dm(x), state) for x in c_mom_eigen[1]])
        for state in psi0:
            probs0_pos.append([expect(ket2dm(x), state) for x in c_pos_eigen[1]])
            probs0_mom.append([expect(ket2dm(x), state) for x in c_mom_eigen[1]])
        for state in psi1:
            probs1_pos.append([expect(ket2dm(x), state) for x in c_pos_eigen[1]])
            probs1_mom.append([expect(ket2dm(x), state) for x in c_mom_eigen[1]])
        for state in psi2:
            probs2_pos.append([expect(ket2dm(x), state) for x in c_pos_eigen[1]])
            probs2_mom.append([expect(ket2dm(x), state) for x in c_mom_eigen[1]])


        # Check that probabilities are normalised (they should sum to 1 or 50 depending on the form)

        diffs_pos = []
        diffs_mom = []
        fisher_info_pos = array([])
        fisher_info_mom = array([])

        print "Calculating differentials"
        # Cycle through the probs and work out the differentials
        # Use a 2nd order central difference theorem
        for timestep in range(0, len(probs0_pos)):
            diffs_pos.append([(-probs2_pos[timestep][xvalue] + 8*probs1_pos[timestep][xvalue] - 8*probsm1_pos[timestep][xvalue] + probsm2_pos[timestep][xvalue])/(12*self.h) for xvalue in range(0, len(probs0_pos[0]))])
            diffs_mom.append([(-probs2_mom[timestep][xvalue] + 8*probs1_mom[timestep][xvalue] - 8*probsm1_mom[timestep][xvalue] + probsm2_mom[timestep][xvalue])/(12*self.h) for xvalue in range(0, len(probs0_mom[0]))])

        # Make into numpy arrays to simplify future operations
        diffs_pos = asarray(diffs_pos, dtype = float64)
        diffs_mom = asarray(diffs_mom, dtype = float64)
        probs0_pos = asarray(probs0_pos, dtype = float64)
        probs0_mom = asarray(probs0_mom, dtype = float64)


        for timestep in range(0, len(diffs_pos)):
            elements_pos = diffs_pos[timestep]*diffs_pos[timestep]/probs0_pos[timestep] # Use array multiplication to calculate each
            elements_mom = diffs_mom[timestep]*diffs_mom[timestep]/probs0_mom[timestep] # Use array multiplication to calculate each
            elements_pos = elements_pos[logical_not(isnan(elements_pos))] #Remove nan from array
            elements_mom = elements_mom[logical_not(isnan(elements_mom))] #Remove nan from array
            fisher_info_pos = append(fisher_info_pos, sum(elements_pos))
            fisher_info_mom = append(fisher_info_mom, sum(elements_mom))

        return [fisher_info_pos, fisher_info_mom]


    def RungeKutta4_Lindblad(self, statevector):
        k1 = self.Lindblad(statevector)
        k2 = self.Lindblad(statevector + (self.h/2.)*k1)
        k3 = self.Lindblad(statevector + (self.h/2.)*k2)
        k4 = self.Lindblad(statevector + self.h*k3)
        return (self.h/6.)*(k1 + 2.*k2 + 2.*k3 + k4)


    def Lindblad(self, statevec): # This function takes a vector consisting of the state and its derivative with respect to gbar.
        return array([- 1j*(self.H.dot(statevec[0]) - statevec[0].dot(self.H)) + self.Lind.dot(statevec[0].dot(self.Linddag)) - (1/2.)*(statevec[0].dot(self.LinddagLind) + self.LinddagLind.dot(statevec[0])),  1j*(self.dH.dot(statevec[0]) - statevec[0].dot(self.dH)) - 1j*(self.H.dot(statevec[1]) - statevec[1].dot(self.H)) + self.Lind.dot(statevec[1].dot(self.Linddag)) - (1/2.)*(statevec[1].dot(self.LinddagLind) + self.LinddagLind.dot(statevec[1]))])


if __name__ == "__main__":
    freeze_support()
    system = Solver('config.yaml')
