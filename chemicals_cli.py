#!/usr/bin/env python3
import sys
import textwrap
import re
from pprint import pprint

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

def is_valid_cas(cas):
    """Check if a string appears to be a CAS number format"""
    pattern = r'^\d{1,7}-\d{2}-\d$'
    return bool(re.match(pattern, cas))

def determine_polarity(cas):
    """Determine if a chemical is polar, non-polar, or in between based on its dipole moment"""
    try:
        dipole_value = dipole.dipole_moment(cas)
        if dipole_value is None:
            return "Unknown (dipole moment data not available)"
        
        # Classify the polarity based on the dipole moment
        if dipole_value < 0.5:
            return "Non-polar (dipole moment < 0.5 debye)"
        elif dipole_value < 1.5:
            return "Slightly polar (dipole moment 0.5-1.5 debye)"
        elif dipole_value < 3.0:
            return "Moderately polar (dipole moment 1.5-3.0 debye)"
        else:
            return "Strongly polar (dipole moment > 3.0 debye)"
    except Exception as e:
        return f"Unknown (error: {e})"

def lookup_chemical(cas_or_name):
    """Look up a chemical by CAS number or name and return all available properties."""
    try:
        # Try to get chemical info
        cas = None
        name = None
        
        # If it's a CAS, use it directly
        if is_valid_cas(cas_or_name):
            try:
                # Verify it exists by looking up molecular weight
                cas = cas_or_name
                MW = identifiers.MW(cas)
                metadata = identifiers.search_chemical(cas)
                if metadata and metadata.common_name:
                    name = metadata.common_name
                else:
                    name = cas  # Use CAS as name if we can't find common name
            except Exception as e:
                print(f"Could not find chemical with CAS: {cas_or_name}, error: {e}")
                return
        else:
            # Try to get the CAS from the name
            try:
                cas = identifiers.CAS_from_any(cas_or_name)
                if cas:
                    name = cas_or_name
                    # Try to get a better name from the database
                    try:
                        metadata = identifiers.search_chemical(cas)
                        if metadata and metadata.common_name:
                            name = metadata.common_name
                    except:
                        pass
            except Exception as e:
                print(f"Could not find chemical: {cas_or_name}, error: {e}")
                return
                
        print(f"\nChemical Information for: {name} (CAS: {cas})")
        print("=" * 50)
        
        results = {}
        
        # Basic identifiers
        try:
            results["identifiers"] = {
                "cas": cas,
                "name": name,
                "formula": None,
                "smiles": None,
                "inchi": None,
                "pubchem_id": None
            }
            
            # Try to get metadata from the search functionality
            try:
                metadata = identifiers.search_chemical(cas)
                if metadata:
                    results["identifiers"]["formula"] = metadata.formula
                    results["identifiers"]["smiles"] = metadata.smiles
                    results["identifiers"]["inchi"] = metadata.InChI
                    results["identifiers"]["pubchem_id"] = metadata.pubchemid
            except:
                pass
        except Exception as e:
            pass
        
        # Critical properties
        try:
            results["critical"] = {
                "Tc": critical.Tc(cas),
                "Pc": critical.Pc(cas),
                "Vc": critical.Vc(cas),
                "Zc": critical.Zc(cas)
            }
        except Exception as e:
            pass
        
        # Thermal properties
        try:
            results["thermal"] = {}
            
            # Individual properties with try/except for each
            try:
                if hasattr(thermal_conductivity, 'Sheffy_Johnson'):
                    results["thermal"]["thermal_conductivity_liquid"] = thermal_conductivity.Sheffy_Johnson(cas, T=298.15)
            except:
                pass
                
            try:
                if hasattr(thermal_conductivity, 'DIPPR9B'):
                    results["thermal"]["thermal_conductivity_gas"] = thermal_conductivity.DIPPR9B(cas, T=298.15)
            except:
                pass
                
            try:
                results["thermal"]["viscosity_liquid"] = viscosity.viscosity_liquid(cas, T=298.15)
            except:
                pass
                
            try:
                results["thermal"]["viscosity_gas"] = viscosity.viscosity_gas(cas, T=298.15)
            except:
                pass
                
            try:
                results["thermal"]["heat_capacity_liquid"] = heat_capacity.Cp_liquid(cas, T=298.15)
            except:
                pass
                
            try:
                results["thermal"]["heat_capacity_gas"] = heat_capacity.Cp_gas(cas, T=298.15)
            except:
                pass
                
            try:
                results["thermal"]["heat_capacity_solid"] = heat_capacity.Cp_solid(cas, T=298.15)
            except:
                pass
                
            # Remove if empty
            if not results["thermal"]:
                del results["thermal"]
                
        except Exception as e:
            pass
        
        # Phase properties
        try:
            results["phase"] = {}
            
            try:
                results["phase"]["Tb"] = phase_change.Tb(cas)
            except:
                pass
                
            try:
                results["phase"]["Tm"] = phase_change.Tm(cas)
            except:
                pass
                
            try:
                results["phase"]["Hfus"] = phase_change.Hfus(cas)
            except:
                pass
                
            try:
                if hasattr(phase_change, 'Hvap'):
                    results["phase"]["Hvap"] = phase_change.Hvap(cas, T=298.15)
            except:
                pass
                
            try:
                results["phase"]["surface_tension"] = interface.sigma(cas, T=298.15)
            except:
                pass
                
            try:
                results["phase"]["vapor_pressure"] = vapor_pressure.Psat(cas, T=298.15)
            except:
                pass
                
            try:
                results["phase"]["triple_temperature"] = triple.Tt(cas)
            except:
                pass
                
            try:
                results["phase"]["triple_pressure"] = triple.Pt(cas)
            except:
                pass
                
            # Remove if empty
            if not results["phase"]:
                del results["phase"]
                
        except Exception as e:
            pass
        
        # Physical properties
        try:
            results["physical"] = {}
            
            try:
                results["physical"]["MW"] = identifiers.MW(cas)
            except:
                pass
                
            try:
                if hasattr(volume, 'density_liquid'):
                    results["physical"]["density_liquid"] = volume.density_liquid(cas, T=298.15)
            except:
                pass
                
            try:
                if hasattr(volume, 'density_solid'):
                    results["physical"]["density_solid"] = volume.density_solid(cas, T=298.15)
            except:
                pass
                
            try:
                if hasattr(volume, 'volume_liquid'):
                    results["physical"]["volume_liquid"] = volume.volume_liquid(cas, T=298.15)
            except:
                pass
                
            try:
                results["physical"]["acentric_factor"] = acentric.omega(cas)
            except:
                pass
                
            # Add dipole moment and polarity classification together
            try:
                dipole_value = dipole.dipole_moment(cas)
                results["physical"]["dipole_moment"] = dipole_value
                if dipole_value is not None:
                    results["physical"]["polarity"] = determine_polarity(cas)
            except:
                pass
                
            try:
                results["physical"]["refractive_index"] = refractivity.refractive_index(cas, T=298.15)
            except:
                pass
                
            try:
                results["physical"]["permittivity"] = permittivity.permittivity(cas, T=298.15)
            except:
                pass
                
            # Remove if empty
            if not results["physical"]:
                del results["physical"]
                
        except Exception as e:
            pass
        
        # Safety and environmental
        try:
            results["safety"] = {}
            
            try:
                results["safety"]["LFL"] = safety.LFL(cas)
            except:
                pass
                
            try:
                results["safety"]["UFL"] = safety.UFL(cas)
            except:
                pass
                
            try:
                results["safety"]["flash_point"] = safety.Tflash(cas)
            except:
                pass
                
            try:
                results["safety"]["autoignition_temp"] = safety.Tautoignition(cas)
            except:
                pass
                
            try:
                results["safety"]["GWP"] = environment.GWP(cas)
            except:
                pass
                
            try:
                results["safety"]["ODP"] = environment.ODP(cas)
            except:
                pass
                
            try:
                results["safety"]["logP"] = environment.logP(cas)
            except:
                pass
                
            # Remove if empty
            if not results["safety"]:
                del results["safety"]
                
        except Exception as e:
            pass
        
        # Check if it's an element in the periodic table - by symbol
        try:
            element_info = None
            if len(cas_or_name) <= 2 and cas_or_name.isalpha():
                # Try to look up by symbol
                try:
                    element_info = elements.periodic_table[cas_or_name]
                except:
                    pass
            
            # Check by name
            if element_info is None:
                try:
                    element_info = elements.periodic_table[name]
                except:
                    pass
                    
            if element_info:
                results["element"] = {
                    "atomic_number": element_info.number,
                    "period": element_info.period,
                    "group": element_info.group,
                    "atomic_weight": element_info.MW
                }
        except Exception as e:
            pass
            
        # Molecular properties
        try:
            results["molecular"] = {}
            
            try:
                if hasattr(miscdata, 'atom_count'):
                    results["molecular"]["atom_count"] = miscdata.atom_count(cas)
            except:
                pass
                
            try:
                results["molecular"]["molecular_diameter"] = lennard_jones.molecular_diameter(cas)
            except:
                pass
                
            try:
                if hasattr(miscdata, 'has_hydroxyl'):
                    results["molecular"]["has_hydroxyl"] = miscdata.has_hydroxyl(cas)
            except:
                pass
                
            try:
                if hasattr(miscdata, 'has_ether'):
                    results["molecular"]["has_ether"] = miscdata.has_ether(cas)
            except:
                pass
                
            # Remove if empty
            if not results["molecular"]:
                del results["molecular"]
        except Exception as e:
            pass
            
        # Print the results
        for category, props in results.items():
            print(f"\n{category.upper()}:")
            for prop, value in props.items():
                if value is not None:
                    # Format value for better readability
                    if isinstance(value, float):
                        formatted_value = f"{value:.6g}"
                    elif isinstance(value, list) and len(value) > 5:
                        formatted_value = f"{value[:5]} ... ({len(value)} items)"
                    else:
                        formatted_value = str(value)
                    
                    # Format the property name for better readability
                    formatted_prop = " ".join(word.capitalize() for word in prop.split("_"))
                    
                    print(f"  {formatted_prop}: {formatted_value}")
                    
        return results
                    
    except Exception as e:
        print(f"Error retrieving chemical data: {e}")
        import traceback
        traceback.print_exc()
        return None

def search_chemicals(query, limit=10):
    """Search for chemicals by partial name using PubChem database"""
    matches = []
    
    # Try to load the database and retrieve chemicals
    db = identifiers.get_pubchem_db()
    
    # Different approaches based on db attributes
    if hasattr(db, 'common_name_index'):
        # Search through common names
        for name, metadata in db.common_name_index.items():
            if query.lower() in name.lower():
                matches.append((name, metadata.CASs))
    
        # If no matches in common names, try synonyms
        if not matches and hasattr(db, 'synonym_index'):
            for synonym, metadata in db.synonym_index.items():
                if query.lower() in synonym.lower():
                    matches.append((metadata.common_name, metadata.CASs))
    else:
        # Alternative approach if common_name_index is not available
        try:
            for i in range(len(db.CAS_index)):
                # Get chemical by index
                metadata = db[i]
                if metadata and metadata.common_name and query.lower() in metadata.common_name.lower():
                    matches.append((metadata.common_name, metadata.CASs))
        except:
            # Fallback - try direct CAS lookup
            try:
                cas = identifiers.CAS_from_any(query)
                metadata = identifiers.search_chemical(cas)
                if metadata:
                    matches.append((metadata.common_name, cas))
            except:
                pass
                
    # Remove duplicates and sort
    unique_matches = {}
    for name, cas in matches:
        if cas not in unique_matches:
            unique_matches[cas] = name
    
    # Convert back to list and sort by name length
    matches = [(name, cas) for cas, name in unique_matches.items()]
    matches.sort(key=lambda x: len(x[0]))
    
    if matches:
        print(f"\nFound {len(matches)} chemicals matching '{query}':")
        for i, (name, cas) in enumerate(matches[:limit], 1):
            print(f"{i}. {name} (CAS: {cas})")
        
        if len(matches) > limit:
            print(f"...and {len(matches) - limit} more.")
            
        return matches[:limit]
    else:
        print(f"No chemicals found matching '{query}'")
        
        # Suggest using CAS lookup if no matches found
        print("Try searching by CAS number instead.")
        return []

def property_calculator_menu():
    """Menu for calculating specific properties for any chemical"""
    properties = {
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
        "12": ("Liquid Density", 
               lambda cas: volume.rhol(cas, T=298.15), 
               "kg/m³"),
        "13": ("Gas Density", 
               lambda cas: volume.rhog(cas, T=298.15, P=101325), 
               "kg/m³"),
        "14": ("Liquid Viscosity", 
               lambda cas: viscosity.mu_l(cas, T=298.15), 
               "Pa·s"),
        "15": ("Gas Viscosity", 
               lambda cas: viscosity.mu_g(cas, T=298.15), 
               "Pa·s"),
        "16": ("Liquid Thermal Conductivity", 
               lambda cas: thermal_conductivity.kl(cas, T=298.15), 
               "W/(m·K)"),
        "17": ("Gas Thermal Conductivity", 
               lambda cas: thermal_conductivity.kg(cas, T=298.15), 
               "W/(m·K)"),
        "18": ("Water Solubility", 
               lambda cas: solubility.solubility_parameter(cas), 
               "J^0.5/m^1.5")
    }
    
    while True:
        print("\nProperty Calculator")
        print("------------------")
        for key, (name, _, unit) in properties.items():
            if unit:
                print(f"{key}. {name} [{unit}]")
            else:
                print(f"{key}. {name}")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nSelect property to calculate: ").lower()
        
        if choice == 'back':
            break
            
        if choice in properties:
            cas_or_name = input("Enter chemical name or CAS number: ")
            
            # Validate and convert to CAS if necessary
            cas = None
            if is_valid_cas(cas_or_name):
                cas = cas_or_name
            else:
                try:
                    cas = identifiers.CAS_from_any(cas_or_name)
                except:
                    print(f"Could not find chemical: {cas_or_name}")
                    continue
                    
            try:
                property_name, calculation_func, unit = properties[choice]
                result = calculation_func(cas)
                
                if choice == "11":  # Special handling for polarity classification
                    print(f"\n{property_name} of {cas_or_name} (CAS: {cas}): {result}")
                else:
                    print(f"\n{property_name} of {cas_or_name} (CAS: {cas}): {result:.6g} {unit}")
                    
                    # Additional information for dipole moment
                    if choice == "8":
                        polarity = determine_polarity(cas)
                        print(f"Polarity Classification: {polarity}")
                        
            except Exception as e:
                print(f"Error calculating {property_name}: {e}")
                
                # Suggest using another function if available
                if choice == "6":  # Surface tension
                    print("Note: Try using option 7 (Acentric Factor) which is more widely available.")
                elif choice == "10":  # Heat of Formation
                    print("Note: Try using option 5 (Heat of Fusion) which is more widely available.")
        else:
            print("Invalid choice, please try again.")

def comparison_menu():
    """Menu for comparing properties of multiple chemicals"""
    while True:
        print("\nChemical Property Comparison")
        print("--------------------------")
        print("1. Compare Critical Temperatures")
        print("2. Compare Boiling Points")
        print("3. Compare Flash Points (Safety)")
        print("4. Compare Acentric Factors")
        print("5. Compare Vapor Pressures")
        print("6. Compare GWP (Global Warming Potential)")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nEnter your choice: ").lower()
        
        if choice == 'back':
            break
            
        if choice in ['1', '2', '3', '4', '5', '6']:
            # Get the property function based on choice
            if choice == '1':
                prop_func = lambda cas: critical.Tc(cas)
                prop_name = "Critical Temperature (K)"
            elif choice == '2':
                prop_func = lambda cas: phase_change.Tb(cas)
                prop_name = "Boiling Point (K)"
            elif choice == '3':
                prop_func = lambda cas: safety.Tflash(cas)
                prop_name = "Flash Point (K)"
            elif choice == '4':
                prop_func = lambda cas: acentric.omega(cas)
                prop_name = "Acentric Factor (dimensionless)"
            elif choice == '5':
                prop_func = lambda cas: vapor_pressure.Psat(cas, T=298.15)
                prop_name = "Vapor Pressure at 298.15 K (Pa)"
            elif choice == '6':
                prop_func = lambda cas: environment.GWP(cas)
                prop_name = "Global Warming Potential"
                
            # Get the chemicals to compare
            chemical_input = input("Enter chemical names or CAS numbers separated by commas: ")
            chemicals = [c.strip() for c in chemical_input.split(',')]
            
            # Get the properties and names
            results = []
            for chem in chemicals:
                cas = None
                if is_valid_cas(chem):
                    cas = chem
                else:
                    try:
                        cas = identifiers.CAS_from_any(chem)
                    except:
                        print(f"Could not find chemical: {chem}")
                        continue
                        
                try:
                    value = prop_func(cas)
                    results.append((chem, cas, value))
                except Exception as e:
                    print(f"Error calculating property for {chem}: {e}")
            
            # Sort the results
            if results:
                results.sort(key=lambda x: x[2] if x[2] is not None else float('-inf'))
                
                # Print the comparison
                print(f"\nComparison of {prop_name}:")
                print("-" * 50)
                for name, cas, value in results:
                    if value is not None:
                        print(f"{name} (CAS: {cas}): {value:.6g}")
                    else:
                        print(f"{name} (CAS: {cas}): Data not available")
            else:
                print("No valid results to compare")
        else:
            print("Invalid choice, please try again.")

def reaction_calculator():
    """Calculate properties for chemical reactions"""
    while True:
        print("\nChemical Reaction Calculator")
        print("--------------------------")
        print("1. Calculate Heat of Formation (gas)")
        print("2. Calculate Heat of Formation (liquid)")
        print("3. Calculate Standard Entropy (gas)")
        print("4. Calculate Gibbs Energy of Formation (gas)")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nEnter your choice: ").lower()
        
        if choice == 'back':
            break
            
        if choice in ['1', '2', '3', '4']:
            # Get chemical
            cas_or_name = input("Enter chemical name or CAS number: ")
            
            # Convert name to CAS if needed
            cas = None
            if is_valid_cas(cas_or_name):
                cas = cas_or_name
            else:
                try:
                    cas = identifiers.CAS_from_any(cas_or_name)
                except:
                    print(f"Could not find chemical: {cas_or_name}")
                    continue
            
            # Temperature for the calculation
            if choice == '1' or choice == '2' or choice == '3' or choice == '4':
                T = float(input("Enter temperature in K (default is 298.15): ") or 298.15)
            
            try:
                # Calculate the property based on choice
                if choice == '1':
                    # Heat of formation in gas phase
                    result = reaction.Hfg(cas, T=T)
                    print(f"\nHeat of Formation (gas) for {cas_or_name} at {T} K: {result:.6g} J/mol")
                elif choice == '2':
                    # Heat of formation in liquid phase
                    result = reaction.Hfl(cas, T=T)
                    print(f"\nHeat of Formation (liquid) for {cas_or_name} at {T} K: {result:.6g} J/mol")
                elif choice == '3':
                    # Standard entropy in gas phase
                    result = reaction.S0g(cas, T=T)
                    print(f"\nStandard Entropy (gas) for {cas_or_name} at {T} K: {result:.6g} J/(mol·K)")
                elif choice == '4':
                    # Gibbs energy of formation
                    if hasattr(reaction, 'Gibbs_formation'):
                        result = reaction.Gibbs_formation(cas, T=T)
                        print(f"\nGibbs Energy of Formation for {cas_or_name} at {T} K: {result:.6g} J/mol")
                    else:
                        # Calculate Gibbs energy from heat of formation and entropy if direct function not available
                        hf = reaction.Hfg(cas, T=T)
                        s0 = reaction.S0g(cas, T=T)
                        result = hf - T*s0
                        print(f"\nGibbs Energy of Formation for {cas_or_name} at {T} K: {result:.6g} J/mol (calculated from Hfg and S0g)")
            except Exception as e:
                print(f"Error calculating property: {e}")
                
                # Show available methods if lookup failed
                try:
                    if choice == '1':
                        methods = reaction.Hfg_methods(cas)
                        print(f"Available methods for Heat of Formation (gas): {methods}")
                    elif choice == '2':
                        methods = reaction.Hfl_methods(cas)
                        print(f"Available methods for Heat of Formation (liquid): {methods}")
                    elif choice == '3':
                        methods = reaction.S0g_methods(cas)
                        print(f"Available methods for Standard Entropy (gas): {methods}")
                except:
                    print("No data available for this chemical in any method.")
        else:
            print("Invalid choice, please try again.")

def element_lookup_menu():
    """Menu for looking up element data from the periodic table"""
    print("\nPeriodic Table Element Lookup")
    print("----------------------------")
    
    while True:
        element = input("\nEnter element symbol, name, or atomic number (or 'back' to return): ")
        if element.lower() == 'back':
            break
            
        try:
            # Try to get element info
            element_info = None
            
            # Try by atomic number
            if element.isdigit():
                atomic_number = int(element)
                for e in elements.periodic_table:
                    if e.number == atomic_number:
                        element_info = e
                        break
            
            # Try by symbol (case insensitive)
            if element_info is None and len(element) <= 2:
                try:
                    element_info = elements.periodic_table[element.capitalize()]
                except:
                    pass
            
            # Try by name (case insensitive)
            if element_info is None:
                for e in elements.periodic_table:
                    if element.lower() == e.name.lower():
                        element_info = e
                        break
                        
            if element_info:
                print(f"\nElement Information: {element_info.name} ({element_info.symbol})")
                print("=" * 50)
                print(f"Atomic Number: {element_info.number}")
                print(f"Atomic Weight: {element_info.MW:.6g} g/mol")
                print(f"Period: {element_info.period}")
                print(f"Group: {element_info.group}")
                
                # Get block information
                block = None
                if element_info.number in elements.s_block:
                    block = "s-block"
                elif element_info.number in elements.d_block:
                    block = "d-block"
                elif element_info.number in elements.p_block:
                    block = "p-block"
                elif element_info.number in elements.f_block:
                    block = "f-block"
                    
                if block:
                    print(f"Block: {block}")
                
                # Print other available properties
                if element_info.phase:
                    phase_names = {'s': 'Solid', 'l': 'Liquid', 'g': 'Gas'}
                    print(f"Phase at STP: {phase_names.get(element_info.phase, element_info.phase)}")
                    
                if element_info.Hf is not None:
                    print(f"Heat of Formation: {element_info.Hf:.6g} J/mol")
                    
                if element_info.S0 is not None:
                    print(f"Standard Entropy: {element_info.S0:.6g} J/(mol·K)")
                    
                if element_info.electronegativity is not None:
                    print(f"Electronegativity: {element_info.electronegativity:.3g}")
                    
                if element_info.covalent_radius is not None:
                    print(f"Covalent Radius: {element_info.covalent_radius:.3g} Å")
                    
                if element_info.vdw_radius is not None:
                    print(f"Van der Waals Radius: {element_info.vdw_radius:.3g} Å")
            else:
                print(f"Could not find element: {element}")
                
        except Exception as e:
            print(f"Error looking up element {element}: {e}")

def main_menu():
    """Main menu for the chemical data CLI"""
    while True:
        print("""
Chemical Data CLI
================

1. Search for a chemical by name
2. Get properties for a specific chemical by name or CAS
3. Calculate specific property for any chemical
4. Compare properties of multiple chemicals
5. Reaction property calculator
6. Look up element data from the periodic table
7. Exit
        """)
        
        choice = input("Enter your choice: ")
        
        if choice == '1':
            query = input("Enter search term: ")
            matches = search_chemicals(query)
            
            if matches:
                # Allow user to select one of the search results
                select = input("\nEnter number to select a chemical (or 'back'): ")
                if select.lower() != 'back' and select.isdigit() and 1 <= int(select) <= len(matches):
                    name, cas = matches[int(select) - 1]
                    lookup_chemical(cas)
                    input("\nPress Enter to continue...")
                    
        elif choice == '2':
            cas_or_name = input("Enter chemical name or CAS number: ")
            lookup_chemical(cas_or_name)
            input("\nPress Enter to continue...")
            
        elif choice == '3':
            property_calculator_menu()
            
        elif choice == '4':
            comparison_menu()
            
        elif choice == '5':
            reaction_calculator()
            
        elif choice == '6':
            element_lookup_menu()
            
        elif choice == '7':
            print("Exiting. Thank you for using Chemical Data CLI!")
            break
            
        else:
            print("Invalid choice, please try again.")

if __name__ == "__main__":
    print("Welcome to the Chemical Data CLI!")
    print("This tool provides access to chemical property data.")
    print("NOTE: For best results, use CAS numbers (e.g., 7732-18-5 for water)")
    main_menu()