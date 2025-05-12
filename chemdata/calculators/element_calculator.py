#!/usr/bin/env python3
"""
Element lookup functionality for periodic table data.
"""
import sys

try:
    from chemicals import elements
except ImportError:
    print("Required libraries not found. Make sure chemicals is installed.")
    print("Try: pip install chemicals thermo fluids ht")
    sys.exit(1)

def find_element(element_identifier):
    """
    Look up an element by symbol, name, or atomic number.
    
    Args:
        element_identifier: Symbol, name, or atomic number of the element
        
    Returns:
        Element object if found, None otherwise
    """
    element_info = None
    
    # Try by atomic number
    if str(element_identifier).isdigit():
        atomic_number = int(element_identifier)
        for e in elements.periodic_table:
            if e.number == atomic_number:
                element_info = e
                break
    
    # Try by symbol (case insensitive)
    if element_info is None and len(str(element_identifier)) <= 2:
        try:
            element_info = elements.periodic_table[str(element_identifier).capitalize()]
        except:
            pass
    
    # Try by name (case insensitive)
    if element_info is None:
        for e in elements.periodic_table:
            if str(element_identifier).lower() == e.name.lower():
                element_info = e
                break
                
    return element_info

def get_element_data(element_identifier):
    """
    Get comprehensive data for an element.
    
    Args:
        element_identifier: Symbol, name, or atomic number of the element
        
    Returns:
        Dictionary of element properties if found, None otherwise
    """
    element_info = find_element(element_identifier)
    
    if not element_info:
        return None
        
    # Build a dictionary of element properties
    data = {
        "name": element_info.name,
        "symbol": element_info.symbol,
        "atomic_number": element_info.number,
        "atomic_weight": element_info.MW,
        "period": element_info.period,
        "group": element_info.group,
    }
    
    # Get block information
    if element_info.number in elements.s_block:
        data["block"] = "s-block"
    elif element_info.number in elements.d_block:
        data["block"] = "d-block"
    elif element_info.number in elements.p_block:
        data["block"] = "p-block"
    elif element_info.number in elements.f_block:
        data["block"] = "f-block"
    
    # Add optional properties if available
    if element_info.phase:
        phase_names = {'s': 'Solid', 'l': 'Liquid', 'g': 'Gas'}
        data["phase_at_stp"] = phase_names.get(element_info.phase, element_info.phase)
        
    if element_info.Hf is not None:
        data["heat_of_formation"] = element_info.Hf
        
    if element_info.S0 is not None:
        data["standard_entropy"] = element_info.S0
        
    if element_info.electronegativity is not None:
        data["electronegativity"] = element_info.electronegativity
        
    if element_info.covalent_radius is not None:
        data["covalent_radius"] = element_info.covalent_radius
        
    if element_info.vdw_radius is not None:
        data["vdw_radius"] = element_info.vdw_radius
    
    return data