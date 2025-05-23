#!/usr/bin/env python3
"""
Property calculator functionality for chemical properties.
"""
import sys

try:
    from chemicals import identifiers, critical, thermal_conductivity, viscosity
    from chemicals import heat_capacity, dipole, interface, phase_change, volume
    from chemicals import triple, vapor_pressure, virial, acentric, lennard_jones
    from chemicals import safety, environment, elements, reaction, refractivity
    from chemicals import miscdata, air, temperature, permittivity, solubility
except ImportError:
    print("Required libraries not found. Make sure chemicals is installed.")
    print("Try: pip install chemicals thermo fluids ht")
    sys.exit(1)

from chemdata.utils.validators import is_valid_cas
from chemdata.core.chemical_data import determine_polarity

def get_property_calculators():
    """Return a dictionary of property calculators with metadata."""
    return {
        "1": ("Critical Temperature", 
              lambda cas: critical.Tc(cas), 
              "K"),
        "2": ("Critical Pressure", 
              lambda cas: critical.Pc(cas), 
              "Pa"),
        "3": ("Boiling Point", 
              lambda cas: phase_change.Tb(cas), 
              "K"),
        "4": ("Melting Point", 
              lambda cas: phase_change.Tm(cas), 
              "K"),
        "5": ("Heat of Fusion", 
              lambda cas: phase_change.Hfus(cas), 
              "J/mol"),
        "6": ("Surface Tension", 
              lambda cas: interface.sigma(cas, T=298.15), 
              "N/m"),
        "7": ("Acentric Factor", 
              lambda cas: acentric.omega(cas), 
              "dimensionless"),
        "8": ("Dipole Moment (Polarity)", 
              lambda cas: dipole.dipole_moment(cas), 
              "debye"),
        "9": ("Molecular Weight", 
              lambda cas: identifiers.MW(cas), 
              "g/mol"),
        "10": ("Heat of Formation (Gas)", 
               lambda cas: reaction.Hfg(cas), 
               "J/mol"),
        "11": ("Polarity Classification",
               lambda cas: determine_polarity(cas),
               ""),
        "12": ("Liquid Volume (COSTALD method)", 
               lambda cas: volume.COSTALD(T=298.15, Tc=critical.Tc(cas), Vc=critical.Vc(cas), omega=acentric.omega(cas)), 
               "m³/mol"),
        "13": ("Liquid Density (Rackett method)", 
               lambda cas: 1.0/(volume.Rackett(T=298.15, Tc=critical.Tc(cas), Pc=critical.Pc(cas), omega=acentric.omega(cas))*identifiers.MW(cas)), 
               "kg/m³"),
        "14": ("Gas Viscosity (Lucas method)", 
               lambda cas: viscosity.Lucas_gas(T=298.15, Tc=critical.Tc(cas), Pc=critical.Pc(cas), Zc=critical.Zc(cas), MW=identifiers.MW(cas), dipole=dipole.dipole_moment(cas)), 
               "Pa·s"),
        "15": ("Liquid Viscosity (Letsou-Stiel method)", 
               lambda cas: viscosity.Letsou_Stiel(T=298.15, MW=identifiers.MW(cas), Tc=critical.Tc(cas), Pc=critical.Pc(cas), omega=acentric.omega(cas)), 
               "Pa·s"),
        "16": ("Liquid Thermal Conductivity (Sheffy-Johnson)", 
               lambda cas: thermal_conductivity.Sheffy_Johnson(T=298.15, MW=identifiers.MW(cas), Tc=critical.Tc(cas), Pc=critical.Pc(cas), omega=acentric.omega(cas)), 
               "W/(m·K)"),
        "17": ("Gas Thermal Conductivity (Eli-Hanley)", 
               lambda cas: thermal_conductivity.Eli_Hanley(T=298.15, MW=identifiers.MW(cas), Tc=critical.Tc(cas), Vc=critical.Vc(cas), Zc=critical.Zc(cas), omega=acentric.omega(cas), dipole=dipole.dipole_moment(cas)), 
               "W/(m·K)"),
        "18": ("Hansen Solubility Parameter (Dispersion)", 
               lambda cas: solubility.hansen_delta_d(cas), 
               "Pa^0.5")
    }

def get_comparison_calculators():
    """Return a dictionary of calculators for property comparison."""
    return {
        '1': {
            'function': lambda cas: critical.Tc(cas),
            'name': "Critical Temperature (K)"
        },
        '2': {
            'function': lambda cas: phase_change.Tb(cas),
            'name': "Boiling Point (K)"
        },
        '3': {
            'function': lambda cas: safety.Tflash(cas),
            'name': "Flash Point (K)"
        },
        '4': {
            'function': lambda cas: acentric.omega(cas),
            'name': "Acentric Factor (dimensionless)"
        },
        '5': {
            'function': lambda cas: vapor_pressure.Psat(cas, T=298.15),
            'name': "Vapor Pressure at 298.15 K (Pa)"
        },
        '6': {
            'function': lambda cas: environment.GWP(cas),
            'name': "Global Warming Potential"
        }
    }

def calculate_property(property_key, cas_or_name):
    """Calculate a specific property for a chemical."""
    properties = get_property_calculators()
    
    if property_key not in properties:
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
        property_name, calculation_func, unit = properties[property_key]
        result = calculation_func(cas)
        return {
            "cas": cas,
            "name": cas_or_name,
            "property_name": property_name, 
            "value": result,
            "unit": unit
        }
    except Exception as e:
        print(f"Error calculating property: {e}")
        
        # Additional information for certain properties
        if property_key == "6":  # Surface tension
            print("Note: Try using acentric factor which is more widely available.")
        elif property_key == "10":  # Heat of Formation
            print("Note: Try using heat of fusion which is more widely available.")
            
        return None