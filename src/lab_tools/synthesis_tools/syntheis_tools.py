# Authorship: {{{
''' 
Dario C. Lewczyk
09_10_25
'''
#}}}
# Imports: {{{
from collections import defaultdict
from typing import Dict, List, Tuple, Callable, Any
import re
from lab_tools.utils import utils
import numpy as np
#}}}
# DopedMaterial: {{{ 
class DopedMaterial:
    # __init__:{{{ 
    def __init__(
            self,
            target_mass = 1,
            x = 0.02,
            precursors:List[Tuple[str, Callable[Any, float], List[str]]] = None,
            printout:bool = True
        ):
        self._target_mass = None
        self._x = None
        self._precursors = None

        self.target_mass = target_mass
        self.x= x
        self.precursors = precursors
        # Calculate the precursor masses
        self.update_precursor_masses(printout=printout)
    #}}}
    # target_mass: {{{
    @property
    def target_mass(self): 
        return self._target_mass
    @target_mass.setter
    def target_mass(self, new_target):
        if not isinstance(new_target, int) and not isinstance(new_target, float):
            raise ValueError('target mass must be a number')
        self._target_mass = new_target 
    #}}} 
    # x: {{{
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, new_x):
        if not isinstance(new_x, float):
            raise ValueError('x must be a float.')
        if new_x < 1 and new_x > 0:
            self._x = new_x 
        else:
            raise ValueError('x must be between 0 and 1')
    #}}}
    # precursors: {{{
    @property
    def precursors(self):
        return self._precursors
    @precursors.setter
    def precursors(self, new_precursors):
        try:
            for formula, func, allowed_elements in new_precursors:
                if not isinstance(formula, str):
                    raise ValueError('first item of tuple must be a str')
                try:
                    func(self.x)
                except:
                    raise ValueError('Second arg in tuple must be a lambda function')
                for el in allowed_elements:
                    if not isinstance(el,  str):
                        raise ValueError('Third argument must be a list of strings')
            self._precursors = new_precursors
        except:
            raise ValueError('Precursor not correctly formatted')
    #}}}
    # update_precursor_masses: {{{
    def update_precursor_masses(self,printout:bool = True,  **kwargs):
        ''' 
        If you want to use this aside from the automatic
        implementation, you can use kwargs.

        target_mass: float or int
        x: float or int 0 < x < 1
        precursors: list of tuples (formula, func, allowed elements)
        '''
        # kwargs: {{{
        self.target_mass = kwargs.get('target_mass', self.target_mass)
        self.x = kwargs.get('x', self.x)
        self.precursors = kwargs.get('precursor', self.precursors)
        #}}}
        self.target_composition, self.precursor_masses = utils.calculate_precursor_masses(
            target_mass= self.target_mass,
            x = self.x,
            precursors = self.precursors
        )
        if printout:
            self.print_precursor_masses()
    #}}}
    # print_precursor_masses: {{{
    def print_precursor_masses(self,):
        ''' 
        This function creates a nice printout to show 
        all of the precursor masses you need to weigh out. 
        '''
        print(f'To get {self.target_mass} g of the composition:\n\t{dict(self.target_composition)} (x = {self.x})')

        for precursor, entry in self.precursor_masses.items(): 
            mass = entry['mass']
            moles = np.format_float_scientific( entry['moles'], 4)
            print(f"\tMass of {precursor} required: {mass:.4f} g ({moles} moles)")
    #}}}
#}}}
