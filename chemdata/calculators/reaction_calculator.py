#!/usr/bin/env python3
"""
Reaction property calculator functionality.
"""
import sys

try:
    from chemicals import identifiers, reaction
except ImportError:
    print("Required libraries not found. Make sure chemicals is installed.")
    print("Try: pip install chemicals thermo fluids ht")
    sys.exit(1)

from chemdata.utils.validators import is_valid_cas

def get_reaction_calculators():
    """Return a dictionary of reaction property calculators with metadata."""
    return {
        "1": {
            "name": "Heat of Formation (gas)",
            "function": lambda cas, T: reaction.Hfg(cas, T=T),
            "unit": "J/mol"
        },
        "2": {
            "name": "Heat of Formation (liquid)",
            "function": lambda cas, T: reaction.Hfl(cas, T=T),
            "unit": "J/mol"
        },
        "3": {
            "name": "Standard Entropy (gas)",
            "function": lambda cas, T: reaction.S0g(cas, T=T),
            "unit": "J/(mol·K)"
        },
        "4": {
            "name": "Gibbs Energy of Formation (gas)",
            "function": lambda cas, T: reaction.Gibbs_formation(cas, T=T) if hasattr(reaction, 'Gibbs_formation') else None,
            "unit": "J/mol"
        }
    }

def calculate_reaction_property(property_key, cas_or_name, temperature=298.15):
    """Calculate a specific reaction property for a chemical."""
    reaction_props = get_reaction_calculators()
    
    if property_key not in reaction_props:
        print(f"Invalid property key: {property_key}")
        return None
    
    # Validate and convert to CAS if necessary
    cas = None
    if is_valid_cas(cas_or_name):
        cas = cas_or_name
    else:
        try:
            cas = identifiers.CAS_from_any(cas_or_name)
        except Exception as e:
            print(f"Could not find chemical: {cas_or_name}, error: {e}")
            return None
    
    try:
        calculator = reaction_props[property_key]
        result = calculator["function"](cas, temperature)
        
        # Special case for Gibbs energy calculation
        if property_key == "4" and result is None:
            # Calculate Gibbs energy from heat of formation and entropy if direct function not available
            hf = reaction.Hfg(cas, T=temperature)
            s0 = reaction.S0g(cas, T=temperature)
            result = hf - temperature*s0
            
        return {
            "cas": cas,
            "name": cas_or_name,
            "property_name": calculator["name"],
            "value": result,
            "unit": calculator["unit"],
            "temperature": temperature
        }
    except Exception as e:
        print(f"Error calculating reaction property: {e}")
        
        # Show available methods if lookup failed
        try:
            if property_key == "1":
                methods = reaction.Hfg_methods(cas)
                print(f"Available methods for Heat of Formation (gas): {methods}")
            elif property_key == "2":
                methods = reaction.Hfl_methods(cas)
                print(f"Available methods for Heat of Formation (liquid): {methods}")
            elif property_key == "3":
                methods = reaction.S0g_methods(cas)
                print(f"Available methods for Standard Entropy (gas): {methods}")
        except:
            print("No data available for this chemical in any method.")
            
        return None