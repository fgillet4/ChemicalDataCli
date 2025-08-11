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

def calculate_vapor_pressure_table(cas, T_start=273.15, T_end=373.15, T_step=10):
    """Calculate vapor pressure table for a chemical over a temperature range using multiple methods."""
    try:
        # Import necessary functions
        from chemicals import vapor_pressure, identifiers, critical, acentric, phase_change
        import numpy as np
        
        results = []
        current_T = T_start
        successful_methods = []
        
        while current_T <= T_end:
            psat = None
            method_used = None
            
            # Method 1: Try the main Psat function (uses best available data)
            try:
                psat = vapor_pressure.Psat(cas, T=current_T)
                if psat is not None and psat > 0:
                    method_used = "Psat (best available)"
            except:
                pass
            
            # Method 2: Try Lee-Kesler correlation if we have critical properties
            if psat is None:
                try:
                    Tc = critical.Tc(cas)
                    Pc = critical.Pc(cas) 
                    omega = acentric.omega(cas)
                    if all(x is not None for x in [Tc, Pc, omega]):
                        psat = vapor_pressure.Lee_Kesler(current_T, Tc, Pc, omega)
                        if psat is not None and psat > 0:
                            method_used = "Lee-Kesler"
                except:
                    pass
            
            # Method 3: Try Ambrose-Walton correlation
            if psat is None:
                try:
                    Tc = critical.Tc(cas)
                    Pc = critical.Pc(cas)
                    omega = acentric.omega(cas)
                    if all(x is not None for x in [Tc, Pc, omega]):
                        psat = vapor_pressure.Ambrose_Walton(current_T, Tc, Pc, omega)
                        if psat is not None and psat > 0:
                            method_used = "Ambrose-Walton"
                except:
                    pass
            
            # Method 4: Try boiling point relation if we have boiling point
            if psat is None:
                try:
                    Tb = phase_change.Tb(cas)
                    Tc = critical.Tc(cas)
                    Pc = critical.Pc(cas)
                    if all(x is not None for x in [Tb, Tc, Pc]):
                        psat = vapor_pressure.boiling_critical_relation(current_T, Tb, Tc, Pc)
                        if psat is not None and psat > 0:
                            method_used = "Boiling-Critical Relation"
                except:
                    pass
            
            # Method 5: Try Sanjari correlation (for refrigerants but may work for others)
            if psat is None:
                try:
                    Tc = critical.Tc(cas)
                    Pc = critical.Pc(cas)
                    omega = acentric.omega(cas)
                    if all(x is not None for x in [Tc, Pc, omega]):
                        psat = vapor_pressure.Sanjari(current_T, Tc, Pc, omega)
                        if psat is not None and psat > 0:
                            method_used = "Sanjari"
                except:
                    pass
            
            # Method 6: Try Edalat correlation
            if psat is None:
                try:
                    Tc = critical.Tc(cas)
                    Pc = critical.Pc(cas)
                    omega = acentric.omega(cas)
                    if all(x is not None for x in [Tc, Pc, omega]):
                        psat = vapor_pressure.Edalat(current_T, Tc, Pc, omega)
                        if psat is not None and psat > 0:
                            method_used = "Edalat"
                except:
                    pass
            
            # Method 7: Try water-specific IAPWS method for water
            if psat is None and cas == "7732-18-5":  # Water CAS number
                try:
                    if 273.15 <= current_T <= 647.096:  # Valid range for IAPWS
                        psat = vapor_pressure.Psat_IAPWS(current_T)
                        if psat is not None and psat > 0:
                            method_used = "IAPWS (water-specific)"
                except:
                    pass
            
            # If we got a valid result, add it to the table
            if psat is not None and psat > 0:
                results.append({
                    'temperature_K': current_T,
                    'temperature_C': current_T - 273.15,
                    'pressure_Pa': psat,
                    'pressure_bar': psat / 100000,
                    'pressure_mmHg': psat * 760 / 101325,
                    'pressure_kPa': psat / 1000,
                    'method': method_used
                })
                
                # Track which methods worked
                if method_used not in successful_methods:
                    successful_methods.append(method_used)
            
            current_T += T_step
        
        # Print information about methods used
        if successful_methods:
            print(f"Successfully calculated vapor pressures using: {', '.join(successful_methods)}")
        
        return results
        
    except Exception as e:
        print(f"Error calculating vapor pressure table: {e}")
        return None

def compare_vapor_pressure_methods(cas, T):
    """Compare all available vapor pressure calculation methods at a single temperature."""
    try:
        from chemicals import vapor_pressure, identifiers, critical, acentric, phase_change
        
        results = {}
        
        # Method 1: Main Psat function
        try:
            psat = vapor_pressure.Psat(cas, T=T)
            if psat is not None and psat > 0:
                results["Psat (best available)"] = psat
        except Exception as e:
            results["Psat (best available)"] = f"Failed: {str(e)}"
        
        # Get critical properties once for reuse
        try:
            Tc = critical.Tc(cas)
            Pc = critical.Pc(cas)
            omega = acentric.omega(cas)
        except:
            Tc = Pc = omega = None
        
        # Method 2: Lee-Kesler
        if all(x is not None for x in [Tc, Pc, omega]):
            try:
                psat = vapor_pressure.Lee_Kesler(T, Tc, Pc, omega)
                if psat is not None and psat > 0:
                    results["Lee-Kesler"] = psat
                else:
                    results["Lee-Kesler"] = "Invalid result"
            except Exception as e:
                results["Lee-Kesler"] = f"Failed: {str(e)}"
        else:
            results["Lee-Kesler"] = "Missing critical properties"
        
        # Method 3: Ambrose-Walton
        if all(x is not None for x in [Tc, Pc, omega]):
            try:
                psat = vapor_pressure.Ambrose_Walton(T, Tc, Pc, omega)
                if psat is not None and psat > 0:
                    results["Ambrose-Walton"] = psat
                else:
                    results["Ambrose-Walton"] = "Invalid result"
            except Exception as e:
                results["Ambrose-Walton"] = f"Failed: {str(e)}"
        else:
            results["Ambrose-Walton"] = "Missing critical properties"
        
        # Method 4: Boiling-Critical Relation
        try:
            Tb = phase_change.Tb(cas)
            if all(x is not None for x in [Tb, Tc, Pc]):
                psat = vapor_pressure.boiling_critical_relation(T, Tb, Tc, Pc)
                if psat is not None and psat > 0:
                    results["Boiling-Critical"] = psat
                else:
                    results["Boiling-Critical"] = "Invalid result"
            else:
                results["Boiling-Critical"] = "Missing boiling point or critical properties"
        except Exception as e:
            results["Boiling-Critical"] = f"Failed: {str(e)}"
        
        # Method 5: Sanjari
        if all(x is not None for x in [Tc, Pc, omega]):
            try:
                psat = vapor_pressure.Sanjari(T, Tc, Pc, omega)
                if psat is not None and psat > 0:
                    results["Sanjari"] = psat
                else:
                    results["Sanjari"] = "Invalid result"
            except Exception as e:
                results["Sanjari"] = f"Failed: {str(e)}"
        else:
            results["Sanjari"] = "Missing critical properties"
        
        # Method 6: Edalat
        if all(x is not None for x in [Tc, Pc, omega]):
            try:
                psat = vapor_pressure.Edalat(T, Tc, Pc, omega)
                if psat is not None and psat > 0:
                    results["Edalat"] = psat
                else:
                    results["Edalat"] = "Invalid result"
            except Exception as e:
                results["Edalat"] = f"Failed: {str(e)}"
        else:
            results["Edalat"] = "Missing critical properties"
        
        # Method 7: IAPWS for water
        if cas == "7732-18-5":  # Water
            try:
                if 273.15 <= T <= 647.096:
                    psat = vapor_pressure.Psat_IAPWS(T)
                    if psat is not None and psat > 0:
                        results["IAPWS (water)"] = psat
                    else:
                        results["IAPWS (water)"] = "Invalid result"
                else:
                    results["IAPWS (water)"] = f"Outside valid range (273.15-647.096 K)"
            except Exception as e:
                results["IAPWS (water)"] = f"Failed: {str(e)}"
        
        return results
        
    except Exception as e:
        print(f"Error comparing vapor pressure methods: {e}")
        return None

def calculate_rotavapor_conditions(cas, T_values, P_values):
    """Calculate vapor pressures for rotavapor conditions (temperature and pressure combinations)."""
    try:
        from chemicals import vapor_pressure, identifiers, critical, acentric, phase_change
        
        results = []
        
        for T in T_values:
            for P in P_values:
                # Calculate vapor pressure at this temperature
                psat = None
                method_used = None
                
                # Try multiple methods (same logic as before but condensed)
                methods_to_try = [
                    ("Psat", lambda: vapor_pressure.Psat(cas, T=T)),
                ]
                
                # Add CSP methods if we have properties
                try:
                    Tc = critical.Tc(cas)
                    Pc = critical.Pc(cas)
                    omega = acentric.omega(cas)
                    if all(x is not None for x in [Tc, Pc, omega]):
                        methods_to_try.extend([
                            ("Lee-Kesler", lambda: vapor_pressure.Lee_Kesler(T, Tc, Pc, omega)),
                            ("Ambrose-Walton", lambda: vapor_pressure.Ambrose_Walton(T, Tc, Pc, omega)),
                            ("Sanjari", lambda: vapor_pressure.Sanjari(T, Tc, Pc, omega)),
                            ("Edalat", lambda: vapor_pressure.Edalat(T, Tc, Pc, omega)),
                        ])
                        
                        # Add boiling point method
                        try:
                            Tb = phase_change.Tb(cas)
                            if Tb is not None:
                                methods_to_try.append(("Boiling-Critical", lambda: vapor_pressure.boiling_critical_relation(T, Tb, Tc, Pc)))
                        except:
                            pass
                except:
                    pass
                
                # Try methods until one works
                for method_name, method_func in methods_to_try:
                    try:
                        psat = method_func()
                        if psat is not None and psat > 0:
                            method_used = method_name
                            break
                    except:
                        continue
                
                if psat is not None:
                    # Determine if evaporation will occur
                    will_evaporate = psat > P
                    evaporation_rate = "High" if psat > P * 2 else "Moderate" if psat > P * 1.2 else "Low" if psat > P else "None"
                    
                    results.append({
                        'temperature_K': T,
                        'temperature_C': T - 273.15,
                        'system_pressure_Pa': P,
                        'system_pressure_mbar': P / 100,
                        'vapor_pressure_Pa': psat,
                        'vapor_pressure_mbar': psat / 100,
                        'pressure_ratio': psat / P,
                        'will_evaporate': will_evaporate,
                        'evaporation_rate': evaporation_rate,
                        'method': method_used
                    })
        
        return results
        
    except Exception as e:
        print(f"Error calculating rotavapor conditions: {e}")
        return None

def create_output_directory(cas):
    """Create organized directory structure for output files."""
    import os
    
    # Create main generated_data directory
    main_dir = "generated_data"
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
        print(f"Created directory: {main_dir}")
    
    # Create CAS-specific subdirectory
    cas_dir = os.path.join(main_dir, cas)
    if not os.path.exists(cas_dir):
        os.makedirs(cas_dir)
        print(f"Created directory: {cas_dir}")
    
    return cas_dir