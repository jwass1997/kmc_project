import fenics as fn
import numpy as np
import warnings
import time
import matplotlib.pyplot as plt

from typing import Union, Dict, OrderedDict

class KMCSimulation():
    def __init__(self,
                 x_scale: Union[int, float],
                 y_scale: Union[int, float],
                 **pars: Dict):

        self.x_scale = x_scale
        self.y_scale = y_scale

        self.num_hops = int(pars["num_hops"])
        self.N_a = pars["num_acceptors"]
        self.N_d = pars["num_donors"]
        self.electrode_config = pars["electrode_config"]
        self.num_electrodes = len(self.electrode_config["electrode_positions"])
        self.num_control_electrodes = len(
            self.electrode_config["control_electrode_positions"])

        self.N = self.N_a + self.num_electrodes + self.num_control_electrodes
        self.electrodes = np.zeros((self.num_electrodes, 3))
        self.control_electrodes = np.zeros((self.num_control_electrodes, 3))
        self.occupations = np.zeros((self.N_a))
        self.occupation_time = np.zeros((self.N_a))
        self.electrode_occupations = np.zeros(
            (self.num_electrodes+self.num_control_electrodes))
        self.hopping_counts = np.zeros((self.N, self.N))
        self.distances = np.zeros((self.N, self.N))
        self.const_site_energies = np.zeros((self.N, ))
        self.site_energies = np.zeros((self.N, ))

        self.resolution = pars["resolution"]
        self.nu0 = pars["nu0"]
        self.a = pars["a"]
        self.kbT = pars["kb"]*pars["T"]
        self.mu = 0
        self.R = (self.N_a / (self.x_scale*self.y_scale))**(-1/2)
        self.A0 = 100*self.kbT

    def initialize_simulation(self):
        ''' Configure electrodes and solve poisson eq. '''

        for i in range(self.num_electrodes):
            self.electrodes[i,
                            0] = self.electrode_config["electrode_positions"][i][0]
            self.electrodes[i,
                            1] = self.electrode_config["electrode_positions"][i][1]
            self.electrodes[i,
                            2] = self.electrode_config["electrode_values"][i]

        for i in range(self.num_control_electrodes):
            self.control_electrodes[i,
                                    0] = self.electrode_config["control_electrode_positions"][i][0]
            self.control_electrodes[i,
                                    1] = self.electrode_config["control_electrode_positions"][i][1]
            self.control_electrodes[i,
                                    2] = self.electrode_config["control_electrode_values"][i]
        # Solve poisson eq.
        self.solve_poisson_eq()
        # Setting positions
        self.set_positions()
        # Randomly placing charges
        random_acceptors = np.random.permutation(self.N_a)[:self.N_a-self.N_d]
        self.occupations[random_acceptors] = 1
        # Calculate distances
        self.calc_dist_matrix()

    def kmc_simulation(self):
        ''' KMC simulation '''

        self.initialize_simulation()
        self.calc_const_site_energy()
        # Dynamical and constant part
        dyn_MA_rates = np.zeros((self.N, self.N))
        const_MA_rates = self.calc_const_MA_rates()
        _time = 0
        for _hop in range(self.num_hops):
            print(_hop)
            for i in range(self.N_a):
                acceptor_interaction = 0
                for j in range(self.N_a):
                    if j is not i:
                        acceptor_interaction += (1 - self.occupations[j]) / self.distances[i, j]
                self.site_energies[i] = self.const_site_energies[i]
                self.site_energies[i] -= self.A0*self.R*acceptor_interaction
            # Dynamical part
            for i in range(self.N):
                for j in range(self.N):
                    # electrode-electrode
                    if i >= self.N_a and j >= self.N_a:
                        dyn_MA_rates[i, j] = 0
                    # electrode-acceptor
                    elif i >= self.N_a and j < self.N_a:
                        if self.occupations[j] == 1:
                            dyn_MA_rates[i, j] = 0
                        else:
                            deltaE = self.site_energies[j] - \
                                self.site_energies[i]
                            if deltaE < 0:
                                dyn_MA_rates[i, j] = self.nu0
                            else:
                                dyn_MA_rates[i, j] = self.nu0 * \
                                    np.exp(-deltaE/(self.kbT))
                    # acceptor-electrode
                    elif i < self.N_a and j >= self.N_a:
                        if self.occupations[i] == 0:
                            dyn_MA_rates[i, j] = 0
                        else:
                            deltaE = self.site_energies[j] - \
                                self.site_energies[i]
                            if deltaE < 0:
                                dyn_MA_rates[i, j] = self.nu0
                            else:
                                dyn_MA_rates[i, j] = self.nu0 * \
                                    np.exp(-deltaE/self.kbT)
                    # acceptor-acceptor
                    elif i < self.N_a and j < self.N_a:
                        if self.occupations[i] > 0:
                            if self.occupations[j] == 0:
                                deltaE = self.site_energies[j] - self.site_energies[i] - self.A0*self.R / \
                                    self.distances[i, j]
                                if deltaE < 0:
                                    dyn_MA_rates[i, j] = self.nu0
                                else:
                                    dyn_MA_rates[i, j] = self.nu0 * \
                                        np.exp(-deltaE/self.kbT)
                            else:
                                dyn_MA_rates[i, j] = 0
            # Rate catalog
            self.MA_rates = const_MA_rates * dyn_MA_rates
            self.cum_rate_catalog = np.cumsum(self.MA_rates.flatten())
            self.cum_rate_catalog = self.cum_rate_catalog / self.cum_rate_catalog[-1]
            # Drawing event from a uniform distribution
            event_idx = np.random.rand()
            for i in range(self.N*self.N):
                if self.cum_rate_catalog[i] >= event_idx:
                    event_idx = i
                    break
            # Transform event_idx back to position in MA_rates
            i = int(event_idx / self.MA_rates.shape[0])
            j = int(event_idx % self.MA_rates.shape[1])
            rate_idx = (int(i), int(j))
            # Update occupation of states
            if i < self.N_a and j < self.N_a:
                self.occupations[i] = 0
                self.occupations[j] = 1
            if i < self.N_a and j >= self.N_a:
                self.occupations[i] = 0
                self.electrode_occupations[j - self.N_a] += 1
            if i >= self.N_a and j < self.N_a:
                self.electrode_occupations[i - self.N_a] -= 1
                self.occupations[j] = 1
            # Sample and increase time
            _xi = np.random.random_sample() + 1e-12
            t_d = -((np.log(_xi)) / self.cum_rate_catalog[-1])
            _time += t_d
            # Update
            for i in range(self.N_a):
                
                self.occupation_time += t_d
            self.hopping_counts[i, j] += 1

        return OrderedDict(
            time=_time, 
            occupations=self.occupations,
            electrode_occupations=self.electrode_occupations,
            occupation_time=self.occupation_time,
            hopping_counts=self.hopping_counts
        )

    @staticmethod
    def fn_onboundary(x, on_boundary):

        return on_boundary

    def solve_poisson_eq(self):
        ''' Calculates the potential self.V '''

        self.fn_electrode_dict = {}
        delta = self.x_scale/10

        for i in range(self.num_electrodes):
            self.fn_electrode_dict[f"e{i}_x"] = self.electrodes[i, 0]
            self.fn_electrode_dict[f"e{i}_y"] = self.electrodes[i, 1]
            self.fn_electrode_dict[f"e{i}"] = self.electrodes[i, 2]

        for i in range(self.num_control_electrodes):
            self.fn_electrode_dict[f"ec{i}_x"] = self.control_electrodes[i, 0]
            self.fn_electrode_dict[f"ec{i}_y"] = self.control_electrodes[i, 1]
            self.fn_electrode_dict[f"ec{i}"] = self.control_electrodes[i, 2]

        self.fn_electrode_BC = ""

        for i in range(self.num_electrodes):
            if self.electrodes[i, 0] == 0 or self.electrodes[i, 0] == self.x_scale:
                self.fn_electrode_BC += (f"x[0] == e{i}_x && "
                                         f"x[1] >= e{i}_y - {delta} && "
                                         f"x[1] <= e{i}_y + {delta} ? e{i} : ")
            else:
                self.fn_electrode_BC += (f"x[1] == e{i}_y && "
                                         f"x[0] >= e{i}_x - {delta} && "
                                         f"x[0] <= e{i}_x + {delta} ? e{i} : ")

        for i in range(self.num_control_electrodes):
            if self.control_electrodes[i, 0] == 0 or self.control_electrodes[i, 0] == self.x_scale:
                self.fn_electrode_BC += (f"x[0] == ec{i}_x && "
                                         f"x[1] >= ec{i}_y - {delta} && "
                                         f"x[1] <= ec{i}_y + {delta} ? ec{i} : ")
            else:
                self.fn_electrode_BC += (f"x[1] == ec{i}_y && "
                                         f"x[0] >= ec{i}_x - {delta} && "
                                         f"x[0] <= ec{i}_x + {delta} ? ec{i} : ")

        self.fn_electrode_BC += f"{self.mu}"
        self.fn_boundary = fn.Expression(self.fn_electrode_BC,
                                         degree=1,
                                         **self.fn_electrode_dict)
        self.fn_mesh = fn.RectangleMesh(fn.Point(0, 0),
                                        fn.Point(self.x_scale, self.y_scale),
                                        self.resolution,
                                        self.resolution)
        self.fn_functionspace = fn.FunctionSpace(self.fn_mesh, "P", 1)
        self.fn_BC = fn.DirichletBC(self.fn_functionspace,
                                    self.fn_boundary,
                                    self.fn_onboundary)
        self.V = fn.TrialFunction(self.fn_functionspace)
        self.fn_V = fn.TestFunction(self.fn_functionspace)
        self.fn_a = fn.dot(fn.grad(self.V), fn.grad(self.fn_V))*fn.dx
        self.fn_f = fn.Constant(0)
        self.fn_L = self.fn_f*self.fn_V*fn.dx
        self.V = fn.Function(self.fn_functionspace)
        fn.solve(self.fn_a == self.fn_L, self.V, self.fn_BC)
    
    def set_positions(self):
        
        self.acceptor_pos = np.random.rand(self.N_a, 2)
        self.donor_pos = np.random.rand(self.N_d, 2)
        self.all_electrode_pos = np.zeros(
            (self.num_electrodes+self.num_control_electrodes, 2))
        for i in range(self.num_electrodes):
            self.all_electrode_pos[i, 0] = self.electrodes[i, 0]
            self.all_electrode_pos[i, 1] = self.electrodes[i, 1]

        for i in range(self.num_control_electrodes):
            self.all_electrode_pos[i+self.num_electrodes,
                                   0] = self.control_electrodes[i, 0]
            self.all_electrode_pos[i+self.num_electrodes,
                                   1] = self.control_electrodes[i, 1]

    def calc_dist_matrix(self):
        ''' Calculating distance-matrix '''

        _pos = np.concatenate(
            [self.acceptor_pos, self.all_electrode_pos], axis=0)
        diffs = _pos[:, np.newaxis, :] - _pos[np.newaxis, :, :]
        self.distances = np.sqrt(np.sum(diffs**2, axis=-1))

    def calc_dist(self, pos_i, pos_j):
        ''' Function to calculate the distance between to positions '''

        distance_ij = np.sqrt(
            (pos_i[0] - pos_j[0])**2 + (pos_i[1] - pos_j[1])**2)

        return distance_ij

    def calc_const_site_energy(self):
        ''' Calculates constant part of the energy of every site '''

        for i in range(self.N_a):
            inv_distances = 1.0 / (np.sqrt((self.donor_pos[:, 0]-self.acceptor_pos[i, 0])**2
                                         + (self.donor_pos[:, 1]-self.acceptor_pos[i, 1])**2))
            self.const_site_energies[i] += self.V(self.acceptor_pos[i, 0], self.acceptor_pos[i, 1])
            self.const_site_energies[i] += self.A0*self.R*np.sum(inv_distances, axis=0)
        self.site_energies[self.N_a: self.N_a +
                           self.num_electrodes] = self.electrodes[:, 2]
        self.site_energies[self.N_a+self.num_electrodes: self.N_a+self.num_electrodes +
                           self.num_control_electrodes] = self.control_electrodes[:, 2]

    def update_sol(self):
        ''' Updating the solution self.V after changing electrode voltages '''

        for i in range(len(self.num_electrodes)):
            self.fn_electrode_dict[f"e{i}"] = self.electrodes[i, 2]
        for i in range(len(self.num_control_electrodes)):
            self.fn_electrode_dict[f"ec{i}"] = self.control_electrodes[i, 2]
        self.fn_boundary = fn.Expression(self.fn_electrode_BC,
                                         degree=1,
                                         **self.fn_electrode_dict)
        self.fn_BC = fn.DirichletBC(self.fn_functionspace,
                                    self.fn_boundary,
                                    self.fn_onboundary)
        fn.solve(self.fn_a == self.fn_L, self.V, self.fn_BC)
        self.calculate_site_energy()

    def calc_const_MA_rates(self):
        ''' Calculates the constant part of the transition rate '''
        
        _const_MA_rate = self.nu0*np.exp(-2*self.distances / (self.a*self.R))
        for i in range(_const_MA_rate.shape[0]):
            for j in range(_const_MA_rate.shape[0]):
                if i == j:
                    _const_MA_rate[i, j] = 0

        return _const_MA_rate

if __name__ == "__main__":

    electrode_config = {
        "electrode_positions": [[0.5, 1.0],],
        "control_electrode_positions": [[0.5, 0.0],],
        "electrode_values": [10],
        "control_electrode_values": [-10]
    }

    pars = {
        "num_hops": 1e4,
        "num_acceptors": 10,
        "num_donors": 3,
        "electrode_config": electrode_config
    }

    # Start time
    start = time.time()
    # Make KMC simulation
    sim_kmc = KMCSimulation(x_scale=1, y_scale=1, **pars)
    results = sim_kmc.kmc_simulation()
    # End time
    end = time.time()
    elapsed_time = end - start
    # Print out real time
    print(f"{elapsed_time} seconds")
    # Total simulation time
    _time=results["time"]
    # Show current
    """ visualize_current(hopping_counts=results["hopping_counts"],
                      acceptor_pos=sim_kmc.acceptor_pos,
                      donor_pos=sim_kmc.donor_pos,
                      total_time=_time
                      ) """
    print(results["hopping_counts"])
    # Plot Potential with Dopants
    x = np.arange(0, 1, 0.01)
    y = np.arange(0, 1, 0.01)
    V_plot = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            V_plot[i, j] = sim_kmc.V(x[i], x[j])

    plot_kwargs = dict(s=50, alpha=.7)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    V_show = ax.imshow(V_plot.transpose(),
                       cmap="seismic",
                       interpolation="bicubic",
                       origin="lower",
                       extent=(0, 1, 0, 1))
    cbar = fig.colorbar(V_show)
    ax.scatter(*sim_kmc.acceptor_pos.T, c="k", **plot_kwargs)
    ax.scatter(*sim_kmc.donor_pos.T, c="orange", **plot_kwargs)
    plt.show()