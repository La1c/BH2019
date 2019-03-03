from isdhic import HiCParser
from isdhic.chromosome import ChromosomeSimulation
import sys
import isdhic
import numpy as np

from scipy import interpolate

from isdhic.core import take_time

from scipy import optimize


"""
This script sets up the posterior probability of the structure of
the X chromosome based on single-cell Hi-C contact data. 
"""

## Each bead represents 50 Kb of chromatin and has a diameter of
## roughly 200 nm (assuming a chromatin density of 12 Mb / mu^3)

resolution  =     500 * 1e3
chrsize     = 22430 * 1e3
n_particles = int(chrsize / resolution)
print n_particles

## Linearly ramped Lennard-Jones potential preventing inter-bead
## clashes. The bead diameter is 4 and chosen for compatibility with
## settings used in protein simulations

forcefield  = 'rosetta'
diameter    = 4.

## Single-cell Hi-C data from Nagano et al., the contact distance
## is factor x diameter which equals to 6. Data are from cell 1

######filename    = '../data/test.txt'
factor      = 1.5

## Beads-on-string-model: consecutive beads are connected by a
## harmonic potential with a flat bottom. The harmonic force
## becomes active when the distance between beads (i, i+1) exceeds
## the bead diameter. The force constant controls the strength of
## that force

k_backbone  = 250.

## Radius of gyration (rog) restraint derived from FISH data. The
## size of the X chromosome is roughly 2.58 x rog. Nagano et al.
## determined a size of 3.7 +/- 0.3 mu based on 10 FISH measurements

n_rog       = 10
mu_rog      = 3.7e3 / (2.58 * (resolution / 1200) / diameter)
sigma_rog   = 0.3e3 / (2.58 * (resolution / 1200) / diameter)
tau_rog     = n_rog / sigma_rog**2

## Read data and map chromosomal positions onto 500Kb beads and
## remove contacts arising from loci that are close in sequence and
## were mapped to the same bead.

parser      = HiCParser(filename, 'X', 'X')
datasets    = parser.parse()
dataset     = datasets[('X','X')]
model       = ('logistic', 'relu')[0]

dataset.coarsen(n_particles, chrsize)

dataset.remove_redundant_contacts()
dataset.remove_self_contacts()

## Set up posterior probability using the above settings

simulation  = ChromosomeSimulation(n_particles,
                                   forcefield = forcefield,
                                   k_backbone = k_backbone,
                                   diameter   = diameter,
                                   factor     = factor,
                                   contact_model = model)

posterior = simulation.create_chromosome(list(dataset))
universe  = simulation.universe
coords    = simulation.params['coordinates']
forces    = simulation.params['forces']

posterior['rog'].data[0] = mu_rog
posterior['rog'].tau     = tau_rog


"""
Inferential structure determination of the X chromosome at 500 kb resolution
using Hamiltonian Monte Carlo.
"""


class HamiltonianMonteCarlo(isdhic.HamiltonianMonteCarlo):

    def next(self):

        result = super(HamiltonianMonteCarlo, self).next()

        if len(self.history) and not len(self.history) % 20:
            print '{0}, stepsize = {1:.3e}, -log_prob = {2:.3e}'.format(
                self.history, self.stepsize, self.state.potential_energy)

        return result

if __name__ == '__main__':

    ## set up X chromosome simulation at 500 kb / 50 kb resolution

    #resolution = 500  
    #filename   = './chrX_cell1_{0}kb.py'.format(resolution)

    #with open(filename) as script:
    #    exec script

    ## start from stretched out chromosome structure

    extended = np.multiply.outer(np.arange(n_particles), np.eye(3)[0]) * diameter
    coords.set(extended)

    ## use Hamiltonian Monte Carlo to sample X chromosome structures from the
    ## posterior distribution

    n_steps  = 1e3                                    ## number of HMC iterations
    n_leaps  = 1e2                                    ## number of leapfrog integration steps
    stepsize = 1e-3                                   ## initial integration stepsize
    
    hmc = HamiltonianMonteCarlo(posterior,stepsize=stepsize)
    hmc.leapfrog.n_steps = int(n_leaps)
    hmc.adapt_until      = int(0.5 * n_steps) * 10
    hmc.activate()

    posterior['contacts'].alpha = 100.
    hmc.stepsize = 1e-3
    
    samples = []

    counter = 0
    with take_time('running HMC'):
        while counter < n_steps:
            samples.append(next(hmc))
            counter += 1


###FINE SMAPLES

"""
This script sets up the posterior probability of the structure of
the X chromosome based on single-cell Hi-C contact data. 
"""

## Each bead represents 50 Kb of chromatin and has a diameter of
## roughly 200 nm (assuming a chromatin density of 12 Mb / mu^3)

resolution  =     50 * 1e3
#chrsize     = 22420 * 1e3
n_particles = int(chrsize / resolution)

## Linearly ramped Lennard-Jones potential preventing inter-bead
## clashes. The bead diameter is 4 and chosen for compatibility with
## settings used in protein simulations

forcefield  = 'rosetta'
diameter    = 4.

## Single-cell Hi-C data from Nagano et al., the contact distance
## is factor x diameter which equals to 6. Data are from cell 1

##########filename    = '../data/test.txt'
factor      = 1.5

## Beads-on-string-model: consecutive beads are connected by a
## harmonic potential with a flat bottom. The harmonic force
## becomes active when the distance between beads (i, i+1) exceeds
## the bead diameter. The force constant controls the strength of
## that force

k_backbone  = 250.

## Radius of gyration (rog) restraint derived from FISH data. The
## size of the X chromosome is roughly 2.58 x rog. Nagano et al.
## determined a size of 3.7 +/- 0.3 mu based on 10 FISH measurements

n_rog       = 10
mu_rog      = 3.7e3 / (2.58 * (resolution / 1200) / diameter)
sigma_rog   = 0.3e3 / (2.58 * (resolution / 1200) / diameter)
tau_rog     = n_rog / sigma_rog**2

## Read data and map chromosomal positions onto 500Kb beads and
## remove contacts arising from loci that are close in sequence and
## were mapped to the same bead.

parser      = HiCParser(filename, 'X', 'X')
datasets    = parser.parse()
dataset     = datasets[('X','X')]
model       = ('logistic', 'relu')[0]

dataset.coarsen(n_particles, chrsize)

dataset.remove_redundant_contacts()
dataset.remove_self_contacts()

## Set up posterior probability using the above settings

simulation  = ChromosomeSimulation(n_particles,
                                   forcefield = forcefield,
                                   k_backbone = k_backbone,
                                   diameter   = diameter,
                                   factor     = factor,
                                   contact_model = model)

posterior = simulation.create_chromosome(list(dataset))
universe  = simulation.universe
coords    = simulation.params['coordinates']
forces    = simulation.params['forces']

posterior['rog'].data[0] = mu_rog
posterior['rog'].tau     = tau_rog

"""
Continue HMC sampling at higher resolution.

This script assumes that 'run_hmc.py' has already been executed in
the *same* python session.
"""

def lift_coords(coarse_coords, n_fine):
    """
    Interpolate 3d fiber with a finer sampling.
    """
    tck, u = interpolate.splprep(np.reshape(coarse_coords,(-1,3)).T)    
    coords = interpolate.splev(np.linspace(0,1,n_fine), tck)

    return np.transpose(coords).flatten()

if __name__ == '__main__':
    
    #resolution = 50
    #filename   = './chrX_cell1_{0}kb.py'.format(resolution)

    #with open(filename) as script:
    #    exec script

    coords.set(lift_coords(samples[-1].positions, n_particles))

    ## use Hamiltonian Monte Carlo to sample X chromosome structures from the
    ## posterior distribution

    n_steps  = 1e3                           ## number of HMC iterations
    n_leaps  = 1e1                           ## number of leapfrog integration steps
    stepsize = 1e-7                          ## initial integration stepsize
    
    hmc_fine = HamiltonianMonteCarlo(posterior,stepsize=stepsize)
    hmc_fine.leapfrog.n_steps = int(n_leaps)
    hmc_fine.adapt_until      = int(1e6) #0.5 * n_steps)
    hmc_fine.activate()

    samples_fine = []

    with take_time('running HMC'):
        while len(samples_fine) < n_steps:
            samples_fine.append(hmc_fine.next())
'''
###EXTRAFINE SAMPLE

"""
This script sets up the posterior probability of the structure of
the X chromosome based on single-cell Hi-C contact data. 
"""

## Each bead represents 50 Kb of chromatin and has a diameter of
## roughly 200 nm (assuming a chromatin density of 12 Mb / mu^3)

resolution  =     10 * 1e3
#chrsize     = 22420 * 1e3
n_particles = int(chrsize / resolution)
print n_particles
## Linearly ramped Lennard-Jones potential preventing inter-bead
## clashes. The bead diameter is 4 and chosen for compatibility with
## settings used in protein simulations

forcefield  = 'rosetta'
diameter    = 4.

## Single-cell Hi-C data from Nagano et al., the contact distance
## is factor x diameter which equals to 6. Data are from cell 1

#########filename    = '../data/test.txt'

## Beads-on-string-model: consecutive beads are connected by a
## harmonic potential with a flat bottom. The harmonic force
## becomes active when the distance between beads (i, i+1) exceeds
## the bead diameter. The force constant controls the strength of
## that force

k_backbone  = 250.

## Radius of gyration (rog) restraint derived from FISH data. The
## size of the X chromosome is roughly 2.58 x rog. Nagano et al.
## determined a size of 3.7 +/- 0.3 mu based on 10 FISH measurements

n_rog       = 10
mu_rog      = 3.7e3 / (2.58 * (resolution / 1200) / diameter)
sigma_rog   = 0.3e3 / (2.58 * (resolution / 1200) / diameter)
tau_rog     = n_rog / sigma_rog**2

## Read data and map chromosomal positions onto 500Kb beads and
## remove contacts arising from loci that are close in sequence and
## were mapped to the same bead.

parser      = HiCParser(filename, 'X', 'X')
datasets    = parser.parse()
dataset     = datasets[('X','X')]
model       = ('logistic', 'relu')[0]

dataset.coarsen(n_particles, chrsize)

dataset.remove_redundant_contacts()
dataset.remove_self_contacts()

## Set up posterior probability using the above settings

simulation  = ChromosomeSimulation(n_particles,
                                   forcefield = forcefield,
                                   k_backbone = k_backbone,
                                   diameter   = diameter,
                                   factor     = factor,
                                   contact_model = model)

posterior = simulation.create_chromosome(list(dataset))
universe  = simulation.universe
coords    = simulation.params['coordinates']
forces    = simulation.params['forces']

posterior['rog'].data[0] = mu_rog
posterior['rog'].tau     = tau_rog

"""
Continue HMC sampling at higher resolution.

This script assumes that 'run_hmc.py' has already been executed in
the *same* python session.
"""


def lift_coords(coarse_coords, n_fine):
    """
    Interpolate 3d fiber with a finer sampling.
    """
    tck, u = interpolate.splprep(np.reshape(coarse_coords,(-1,3)).T)    
    coords = interpolate.splev(np.linspace(0,1,n_fine), tck)

    return np.transpose(coords).flatten()

if __name__ == '__main__':
    
    #resolution = 50
    #filename   = './chrX_cell1_{0}kb.py'.format(resolution)

    #with open(filename) as script:
    #    exec script

    coords.set(lift_coords(samples_fine[-1].positions, n_particles))

    ## use Hamiltonian Monte Carlo to sample X chromosome structures from the
    ## posterior distribution

    n_steps  = 1e3                           ## number of HMC iterations
    n_leaps  = 1e1                           ## number of leapfrog integration steps
    stepsize = 1e-7                          ## initial integration stepsize
    
    hmc_superfine = HamiltonianMonteCarlo(posterior,stepsize=stepsize)
    hmc_superfine.leapfrog.n_steps = int(n_leaps)
    hmc_superfine.adapt_until      = int(1e6) #0.5 * n_steps)
    hmc_superfine.activate()

    samples_superfine = []

    with take_time('running HMC'):
        while len(samples_superfine) < n_steps:
            samples_superfine.append(hmc_superfine.next())
'''

