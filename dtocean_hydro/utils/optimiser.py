# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#    Copyright (C) 2017-2018 Mathew Topper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
This module contains the classes used to identify the array layout that
maximise the annual energy production of the array, under the specified
constraints.

.. module:: optimiser
   :platform: Windows
   :synopsis: Optimisation algorithms

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import pickle
from math import sqrt

import cma
import numpy as np


class SearchOptimum(object):
    """
    SearchOptimum: the class is used to optimise the array layout for both
        tidal and wave cases

    Args:
        hyd_obj (Tidal or Wave Object): array_hydrodynamic solver object that
            contains the methods to assess the energy production of the array.
        array_obj (Array_pkg): array_pkg object used to generate new layouts
        val (str/numpy.ndarray)[-/m]: depending on the opt input the variable
            will contains a string to identify the parameterised array layout
            or a collection of point to be optimised using
            contraction/expansion
        opt (int)[-]: optimisation option method
        max_num_dev (float) [-]: max number of devices allowed in the array
            estimated from the array rated power and the machine rated power
        min_q_factor (float)[-]: defines the minimum allowed q-factor
        nogo_areas (list) [m]: list containing the UTM coordinates of the nogo
            areas poligons expressed as [X(Northing),Y(Easting)].

        debug (boolean): if set to True, plots and additional command line
            outputs are issued.

    Attributes:
        main (WP2 class): same as args
        problemOrder (int): dimension of the problem to be solved, parameters
            space order
        combineparam (numpy.ndarray): describes the correlation between
            parameters if any. This is used to reduce the parameters space
            order.
        offsetparam (numpy.ndarray): describes the parameters offset if any.
            Used alog with the combineparam attribute
        _debug (boolean): debug flag
        _Val (str/numpy.ndarray): array structure description
        _Opt (int): array type option.
    """
    
    def __init__(self, optim_func,
                       hyd_obj,
                       array_obj,
                       val,
                       opt,
                       max_num_dev,
                       min_q_factor,
                       nogo_areas,
                       debug=False):
        
        self.nogo_areas = nogo_areas
        self._optim_func = optim_func
        self._array = array_obj
        self._hyd_obj = hyd_obj

        self._debug = debug
        self._Val = val
        self._Opt = opt
        self._normalisation_point = max(array_obj.Dmin)

        # set search space bounds
        self._min_bound = 0.
        self._max_bound = 50.

        # set result constraints
        self._min_q_factor = min_q_factor
        self._max_num_dev = max_num_dev
        self._min_dist = min(array_obj.Dmin)

        self.__set_problem_parameters()


    def __set_problem_parameters(self):
        """
        __set_problem_parameters: method used to set up the problem order based
            on the __Opt and _Val attributes
            
        Returns:

        """
        self.opt_dim = 4  # default case
        self.par_a = np.eye(self.opt_dim)
        self.par_b = np.zeros(self.opt_dim)
        
        if self._Opt == 1:
            if self._Val == 'rectangular':
                self.opt_dim = 2
                self.par_a = self.par_a[:, :2]
                self.par_b[-1] = -np.pi/2
            elif self._Val == 'staggered':
                self.opt_dim = 2
                self.par_a = self.par_a[:, :2]*0.
                self.par_a[0,0] = 1.
                self.par_a[1,0] = 1.
                self.par_a[2,1] = 1.
                self.par_a[3,1] = -1.
            elif self._Val == 'full':
                self.par_a[-1,-1] = -1.
            else:
                err_str = ('Error[UserArray input]: For the Option 1 the '
                           'entered Value is incorrect')
                raise ValueError(err_str)
        elif self._Opt == 3:
            self.opt_dim = 2
            self.par_a = self.par_a[:, :2]
        else:
            raise IOError("The specified array option is out of bounds")

    def eval_optimal_layout(self):
        """
        eval_optimal_layout: calls the method specified in the opt_method args
            and send the results back to the calling environment

        Args:
            opt_method_id (int): defines the optimisation method used in the
                                 call.
                                    1 - Evolutionary Strategy
                                    2 - Monte Carlo
                                    3 - Brutal Force

        Returns:
            0 if the optimisation was successful
            -1 if no optimisation was achieved
        """
        # call the optimisation method
        xopt = self._optim_func(self)
        
        # rescale the optimal solution, if any. 
        if xopt == -1:
            errStr=("Error[OptimisationResults]: No array configuration "
                    "satisfies the given optimisation constrains")
            raise ValueError(errStr)
        
        xmap = self.scale_param(xopt)
        
        module_logger.info("Optimal configuration parameters:")
        module_logger.info("Inter-column distance: {}".format(xmap[0]))
        module_logger.info("Inter-row distance: {}".format(xmap[1]))
        module_logger.info("Column angle: {}".format(xmap[2]))
        module_logger.info("Row angle: {}".format(xmap[3]))
        
        if self._Opt == 1:
            NR, NC, IR, IC, beta, psi = self.param_conditioning(xmap)
            self._array.generator(NR, NC, IR, IC, beta, psi)
        else:
            self._array.coord = (self._Val*int(xmap[0]) / 
                                                 self._normalisation_point)
            
        module_logger.info('Ending the optimisation loop.....')
        module_logger.info('Re-generating the best array layout')
    
    # TODO: Staggered doesn't use column spacing?
    def estimate_start_point(self):
        
        def remap(val, sc, IR, IC, beta, psi):
            if val == 'rectangular':
                X_norm = ((float(IR)-sc)/sc, (float(IC)-sc)/sc)
            elif val == 'staggered':
                X_norm = ((float(IR)-sc)/sc, beta/(0.1*(np.pi/2)))
            elif val == 'full':
                X_norm = ((float(IR)-sc)/sc,
                          (float(IR)-sc)/sc,
                          beta/(0.1*(np.pi/2)),
                          beta/(0.1*(np.pi/2)))
            else:
                X_norm = (IR/sc, IC/sc)
            
            if abs(X_norm[0]-5) > 1:
                sc = float(IR)/6
                return remap(val, sc, IR, IC, beta, psi)
            
            return (X_norm, sc)
        
        module_logger.info("Estimating the optimisation starting point")
        max_eval = 100
        
        # arrange a fake staggered grid for full array case
        if self.opt_dim > 2:  
            par_a = np.eye(4)
            par_b = np.zeros(4)
            par_a = par_a[:, :2]*0.
            par_a[0,0] = 1.
            par_a[1,0] = 1.
            par_a[2,1] = 1.
            par_a[3,1] = -1.
        else:
            par_a = self.par_a
            par_b = self.par_b
        
        # un-normalize optimization variable
        x_b = np.array([5,5],'f')
        xmap = np.dot(par_a, x_b)
        xmap[0] *= self._normalisation_point
        xmap[1] *= self._normalisation_point
        xmap[0] += self._normalisation_point
        xmap[1] += self._normalisation_point
        xmap[2] *= 0.1 * np.pi/2
        xmap[3] *= 0.1 * np.pi/2
        
        xNorm = xmap + par_b
        
        NR, NC, IR, IC, beta, psi = self.param_conditioning(xNorm)
        IC = self._array.Dmin[0]+1.01
        IR = self._array.Dmin[1]+1.01
        
        ind = 0
        
        while True:
            
            if self._Opt == 3:
                scaleX = IC/self._normalisation_point
                scaleY = IR/self._normalisation_point
                self._array.coord = self._Val*np.array([scaleX, scaleY])
                self._array.checkMinDist()
            else:
                self._array.generator(NR, NC, IR, IC, beta, psi)
            
            inside = self._array.checkout()
            
            # check conditions prior to solve the array interaction
            if inside.any() and not self._array.minDist_constraint:
                
                if self._array.coord[inside].shape[0] > self._max_num_dev:
                    IC *= 1.05
                    IR *= 1.05
                else:
                    return remap(self._Val,
                                 self._normalisation_point,
                                 IC,
                                 IR, 
                                 beta,
                                 psi)
            
            else:
                
                IC /= 1.05
                IR /= 1.05
            
            ind += 1
            
            if ind > max_eval:
                return (), self._normalisation_point 
        
        raise RuntimeError('Insufficient facts always invite danger.')
    
    def param_conditioning(self, x_scale):
        """
        param_conditioning: the function applies truncation to the scale
            parameters.
        Args:
            x_scale (numpy.ndarray): array of scaled parameters
        
        Returns:
            array_vals (list): conditioned scaled parameters
        """
        Nb = self._max_num_dev
        
        IR = int(x_scale[1])
        IC = int(x_scale[0])
        beta = int(x_scale[2]*100)/100.0
        psi = int(x_scale[3]*100)/100.0
        
        combos = np.array([[0, 0],
                           [0, 1],
                           [1, 0],
                           [1, 1]])
        devs = np.zeros(4)
        
        for i, combo in enumerate(combos):
            
            # the factor 10 is used because for small skewing angles there is a
            # risk to do not fill the lease area. The extra number of bodies 
            # do not affect the calculation.
            NR = 10 * int(sqrt(Nb)) + combo[0]
            NC = 10 * int(sqrt(Nb)) + combo[1]
            
            self._array.generator(NR, NC, IR, IC, beta, psi)
            inside = self._array.checkout(nogo_list=self.nogo_areas)
            n_devs = self._array.coord[inside].shape[0]
            
            devs[i] = n_devs
        
        most_devs_idx = np.argmax(devs)
        bost_combo = combos[most_devs_idx]
        
        NR = 10 * int(sqrt(Nb)) + bost_combo[0]
        NC = 10 * int(sqrt(Nb)) + bost_combo[1]
        
        return NR, NC, IR, IC, beta, psi
    
    def scale_param(self, x_norm):
        """
        scale_param: the function denormalise the parameters.
        Args:
            x_norm (numpy.ndarray): array of normalised parameters

        Returns:
            x_scale (numpy.ndarray): array of scaled parameters
        """
        xmap = np.dot(self.par_a, x_norm)
        xmap[0] *= self._normalisation_point
        xmap[1] *= self._normalisation_point
        xmap[0] += self._normalisation_point
        xmap[1] += self._normalisation_point
        xmap[2] *= 0.1 * np.pi/2
        xmap[3] *= 0.1 * np.pi/2
        x_scale = xmap + self.par_b
        
        return x_scale
         
    def optimCostFunNorm(self, x):
        """
        optimCostFunNorm: the function scale the parameters and calls the
                            optimCostFun method.
        Args:
            x (list): list of normalised parameters

        Returns:
            fval[0] (float): annual energy production for the given array
            fval[1] (float): q factor for the given array
        """
        x_b = np.array(x)
        
        # un-normalize optimization variable
        xNorm = self.scale_param(x_b)

        # call dernorm cost fun
        fval = self.optimCostFun(xNorm)
        return fval[0], fval[1]
    
    def optimCostFun(self, x):
        """
        optimCostFun: the method calculate the AEP and q-factor for the given
            configuration, calling either the tidal or the wave modules
        
        Args:
            x (list): list of 4 parameters used to build the array layout
        
        Return:
            AEP (float): annual energy production for the given array
            q (float): q factor for the given array
        """
        
        NR, NC, IR, IC, beta, psi = self.param_conditioning(x)
        
        if beta < 0.05 and self._Opt != 3:
            module_logger.debug("The angle between rows and "
                                "columns is too close to zero.")
        
        if self._Opt == 3:
            
            scaleX = IC/self._normalisation_point  
            scaleY = IR/self._normalisation_point
            self._array.coord = self._Val*np.array([scaleX, scaleY])
            self._array.checkMinDist()
            
        else:
            
            module_logger.info("Layout parameters: IC: {} IR: {} beta: {} "
                               "psi: {}".format(IC, IR, beta, psi))
            
            self._array.generator(NR, NC, IR, IC, beta, psi)
        
        inside = self._array.checkout(nogo_list=self.nogo_areas)
        if self._debug: self._array.show(inside)
        
        # check conditions prior to solve the array interaction
        if inside.any() and not self._array.minDist_constraint:
            
            n_devs = self._array.coord[inside].shape[0]
            
            if n_devs > self._max_num_dev:
                
                msg_str = ("Layout not valid: {} devices generated, "
                           "but maximum is {}").format(n_devs,
                                                       self._max_num_dev)
                module_logger.info(msg_str)
                
                # return the squared error from actual value and bound
                # due to the computational constraints this penality term is
                # magnified by a factor 2
                maxdev_error = -((self._array.coord[inside].shape[0] / 
                                      self._max_num_dev - 1) * 200) ** 2
                
                return maxdev_error, -1
            
            if self._debug: module_logger.info("OK constr")
            
            # solve the array interaction
            res = self._hyd_obj.energy(self._array.coord[inside])
            
            module_logger.info("Number of devices: {} AEP: {} q-factor: "
                               "{}".format(self._array.coord[inside].shape[0],
                                           res.AEP_array,
                                           res.q_array))
            
            if res.q_array >= self._min_q_factor:
                
                if self._debug:
                    module_logger.info("Valid config: actual q-factor {} --> "
                                        "min q-factor {}".format(
                                                        res.q_array,
                                                        self._min_q_factor))
                
                return res.AEP_array, res.q_array
            
            if self._debug: module_logger.info("Not valid: q < q_min")
            
            # return the squared error from actual value and bound
            qfactor_error = -((res.q_array /
                               self._min_q_factor - 1) * 100.) ** 2
            
            return qfactor_error, res.q_array
        
        if self._debug:
            module_logger.warning('WARNING! For the given configuration '
                                  'no device is inside the active area!. '
                                  'No calculation is performed!')
        
        mindist_error = -((self._array._actual_mindist / 
                                       self._min_dist - 1) * 100.) ** 2
        
        return mindist_error, -1


def method_brutal_force(searcher, N=5):
    
    x = np.linspace(searcher._min_bound,searcher._max_bound,N)
    y = np.linspace(searcher._min_bound,searcher._max_bound,N)
    fit = np.zeros((N*N))
    ind = -1
    for ii, inter_col in enumerate(x):
        for jj, beta in enumerate(y):
            ind += 1
            print('iteration {} over {}'.format(ind+1, len(x)*len(y)))
            fit[ind] = -searcher.optimCostFunNorm((inter_col, beta))[0]
    # index = np.unravel_index(fit.argmin(), fit.shape)
    index = fit.argmin()
    pickle.dump([fit, x, y],
                open("optimisation_results_brutal_force.pkl", "wb"))

    if fit[index] > 0:
        return -1
    else:
        return (x[index//N], y[int(index%N)])


def method_cma_es(searcher,
                  tolfun=1e1,
                  tolx=1e-3,
                  maxiter=200,
                  maxfevals=2000):
    """
    calls the cma package to optimise the power production
        of the array

    Args (optional):
        tolfun (float)[W]: minimun allowed variation of the fit to decide
            for the solution stagnation
        tolx (float)[-]: minimun allowed variation of the parameters to
            decide for the solution stagnation
        maxiter (int)[-]: max number of population regeneration
        maxfevals (int)[-]: max number of total function evaluation

    Returns:
         x (list): list of normalised parameters that represent the best
           solution found
    """
    
    x0, searcher._normalisation_point = searcher.estimate_start_point()
    if not x0:
        warning_str = ('Could not find a suitable starting point '
                       'for the optimiser, the centroid of the '
                       'parameter space is used instead')
        module_logger.warning(warning_str)
        x0 = searcher.opt_dim * [(searcher._min_bound
                                          + searcher._max_bound) / 2.]
    
    es = cma.CMAEvolutionStrategy(
                    x0,
                    2,
                    {'bounds': [searcher._min_bound, searcher._max_bound],
                     'verb_disp': 0})
    
    es.opts.set('tolfun', tolfun)
    es.opts.set('tolx', tolx)
    es.opts.set('maxiter', maxiter)
    es.opts.set('maxfevals',maxfevals)

    while not es.stop():
        
        solutions = es.ask()
        
        # reduce the significant digits of the search space
        # solutions = [np.around(s, decimals=1) for s in solutions]
        temp = [searcher.optimCostFunNorm(s) for s in solutions]
        fitness = [(-el[0]) for el in temp]
        es.tell(solutions, fitness)
        
        if searcher._debug:
            es.logger.add()
            es.disp(10)
    
    if searcher._debug:
        es.result_pretty()
    
    if es.best.f > 0.:
        return -1
    else:
        return (es.best.x).tolist()


def method_monte_carlo(searcher, maxiter=5):
    """
    optimise the array layout using the Motecarlo simulation approach

    Args:
        n_MC (int): number of simulation to be run. Since there is no
            rational behind this method, but everything is based on the
            randomness of the solution, the number of n_MC is directly
            affecting the stability of the solution

    Returns:
        x (list): list of normalised parameters that represent the best
          solution found
    """
    xx = np.random.random_sample((searcher.opt_dim, maxiter))
    xx *= (searcher._max_bound - searcher._min_bound)
    xx += searcher._min_bound

    fit = np.zeros(maxiter)

    for i in range(maxiter):
        module_logger.info('iteration #: {}'.format(i))

        x = xx[:,i]
        fit[i] = -searcher.optimCostFunNorm(x)[0]
    # find max average energy for arrays with q-factor larger than q_min
    index = fit.argmin()
    pickle.dump([fit, xx],
                open("optimisation_results_brutal_force.pkl", "wb"))

    if searcher._debug:
            module_logger.info('AEP for the different configurations:')
            module_logger.info(fit)
            module_logger.info('Optimal configuration features:')
            module_logger.info(xx[:,index].tolist())

    if fit[index] > 0:
        return -1
    else:
        return xx[:,index].tolist()
