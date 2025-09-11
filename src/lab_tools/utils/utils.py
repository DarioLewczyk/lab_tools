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
#}}}
# atomic weights: {{{

atomic_weights = { 
"H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, 
"C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, 
"Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, 
"S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, 
"Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, 
"Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922, 
"Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, 
"Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98, "Ru": 101.07, "Rh": 102.91, 
"Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76, 
"Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91, 
"Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.96, 
"Gd": 157.25, "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, 
"Yb": 173.05, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, 
"Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38, 
"Pb": 207.2, "Bi": 208.98, "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226, 
"Ac": 227, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237, "Pu": 244, "Am": 243, 
"Cm": 247, "Bk": 247, "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259, "Lr": 266, 
"Rf": 267, "Db": 268, "Sg": 269, "Bh": 270, "Hs": 277, "Mt": 278, "Ds": 281, "Rg": 282, 
"Cn": 285, "Nh": 286, "Fl": 289, "Mc": 290, "Lv": 293, "Ts": 294, "Og": 294
}
    
#}}}
# parse_formula: {{{
def parse_formula(formula): 
    # multiply_counts: {{{
    def multiply_counts(counts, factor): 
        return {element: count * factor for element, count in counts.items()} 
    #}}}
    # parse_segment: {{{
    def parse_segment(segment): 
        stack = [] 
        current = defaultdict(int) 
        i = 0 
        # Loop through the formula string: {{{
        while i < len(segment): 
            # Handle cases where a sub formula is in (): {{{
            if segment[i] == '(': 
                stack.append(current) 
                current = defaultdict(int) 
                i += 1 
            #}}}
            # Hanlde the end case for a sub formula in (): {{{
            elif segment[i] == ')': 
                i += 1 
                m = re.match(r'\d*', segment[i:]) 
                multiplier = int(m.group()) if m.group() else 1 
                i += len(m.group()) 
                temp = multiply_counts(current, multiplier) 
                current = stack.pop() 
                for elem, cnt in temp.items(): 
                    current[elem] += cnt 
            #}}}
            # Handle regular symbols: {{{
            else: 
                m = re.match(r'([A-Z][a-z]?)(\d*\.?\d*)', segment[i:])
                if m: 
                    elem = m.group(1) 
                    try:
                        cnt = int(m.group(2)) if m.group(2) else 1 
                    except:
                        cnt = float(m.group(2)) if m.group(2) else 1
                    current[elem] += cnt 
                    i += len(m.group(0)) 
                else: 
                    i += 1 # Skip invalid characters 
            #}}}
        #}}}
        return current 
    #}}}
    # Handle hydrate cases: {{{
    parts = formula.split('*')
    total = defaultdict(int)

    for part in parts:
        # Check for leading multiplier (e.g., 5H2O)
        m = re.match(r'^(\d+)([A-Z].*)$', part)
        if m:
            multiplier = int(m.group(1))
            segment = m.group(2)
            parsed = parse_segment(segment)
            parsed = multiply_counts(parsed, multiplier)
        else:
            parsed = parse_segment(part)
        for elem, cnt in parsed.items():
            total[elem] += cnt
    #}}}
    return dict(total)
#}}}
# molar_mass: {{{
def molar_mass(composition: Dict[str, int]) -> float: 
    """
    Calculates the molar mass of a compound given its elemental composition.
    """ 
    return sum(atomic_weights[element] * count for element, count in composition.items()) 
#}}}
# combine_compositions: {{{
def combine_compositions(weighted_compositions, allowed_elements=None): 
    """ 
    This allows you to make compositions using tuples of substances
    and their molar fractions

    Parameters: 
        - weighted_compositions: list of tuples (composition_dict, weight) 
        - allowed_elements: set of elements to include (optional) 
    Returns: 
        - defaultdict(float) with combined composition 
        """ 
    result = defaultdict(float) 
    for comp_dict, weight in weighted_compositions: 
        for elem, count in comp_dict.items(): 
            if allowed_elements is None or elem in allowed_elements: 
                result[elem] += count * weight 
    return result
#}}}
# calculate_net_charge: {{{
def net_charge(
        composition: Dict[str, float]
        oxidation_states: Dict[str, float],
        polyatomic_ions = { 
            "NO3": {"composition": {"N": 1, "O": 3}, "charge": -1}, 
            "SO4": {"composition": {"S": 1, "O": 4}, "charge": -2}, 
            "PO4": {"composition": {"P": 1, "O": 4}, "charge": -3}, 
            # Add more as needed 
        }
        ) -> float: 
    ''' 
    composition: The comp dictionary

    oxidation_states: a dict of the atoms included in your formula
        and their oxidation states
    polyatomic_ions: a dict of the polyatomic ions and their charge
    '''
    total_charge = 0.0 
    remaining = composition.copy() 
    # Check for polyatomic ions 
    for ion, data in polyatomic_ions.items(): 
        ion_comp = data["composition"] 
        min_ratio = min(remaining.get(elem, 0) / count for elem, count in ion_comp.items()) 
        if min_ratio >= 1: 
            ion_count = int(min_ratio) 
            total_charge += ion_count * data["charge"] 
            for elem, count in ion_comp.items(): 
                remaining[elem] -= ion_count * count 
    # Add remaining atomic charges 
    for elem, count in remaining.items(): 
        total_charge += oxidation_states.get(elem, 0) * count 
    return total_charge
#}}}
# calculate_precursor_masses: {{{
def calculate_precursor_masses(
        target_mass: float, 
        x:float, 
        precursors: List[Tuple[str, Callable[Any, float], List[str]]],
        debug = False
    ) -> Dict[str, float]: 
    """ 
    Calculates the mass of each precursor needed to 
    synthesize the target compound. precursors: 
    List of tuples (precursor_formula, molar_fraction, included_elements) 
    """ 
    # Determine the target composition: {{{
    allowed_elements = set()
    weighted_compositions = []
    for formula, func, elements in precursors:
        if debug:
            print(f'Formula: {formula}\nfunc: {func}\nElements: {elements}')
        # Loop through the precursors and collect unique elements
        allowed_elements.update(elements) # Updates the allowed elements list without duplicates
        comp = parse_formula(formula) # Get the composition dictionary for the formula
        mf = func(x) # Calculate the molar fraction for the comp
        weighted_compositions.append((comp, mf)) # Add the composition
    allowed_elements = list(allowed_elements)
 
    target_composition = combine_compositions(
            weighted_compositions= weighted_compositions,
            allowed_elements=allowed_elements
    )
    if debug: 
        print(f'Target comp: \n\t{target_composition}\n\tAllowed elements: {allowed_elements}')
    #}}}
    # Calculate total molar mass of the target compound 
    target_molar_mass = molar_mass(target_composition) 
    target_moles = target_mass / target_molar_mass 
    precursor_masses = {} 
    for precursor_formula, mf_func, included_elements in precursors: 
        # Calculate the molar fraction: {{{
        try: 
            molar_fraction = mf_func(x) # This 
        except:
            raise ValueError('Must input a lambda function for the second parameter in your precursor tuples')
        #}}}
        precursor_comp = parse_formula(precursor_formula) 
        filtered_comp = {el: count for el, count in precursor_comp.items() if el in included_elements} 
        precursor_molar_mass = molar_mass(filtered_comp) 
        precursor_moles = target_moles * molar_fraction 
        precursor_mass = precursor_moles * precursor_molar_mass 
        precursor_masses[precursor_formula] = {
                'mass': precursor_mass, 
                'moles': precursor_moles
        }
    return target_composition, precursor_masses 
#}}}

