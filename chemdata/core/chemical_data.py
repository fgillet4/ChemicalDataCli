#!/usr/bin/env python3
"""
Core functionality for retrieving and processing chemical data.
"""
import sys
import traceback
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

from chemdata.utils.validators import is_valid_cas

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

def get_all_properties(cas_or_name):
    """Calculate and display all available properties for a chemical"""
    
    # Validate and convert to CAS if necessary
    cas = None
    if is_valid_cas(cas_or_name):
        cas = cas_or_name
    else:
        try:
            cas = identifiers.CAS_from_any(cas_or_name)
        except:
            print(f"Could not find chemical: {cas_or_name}")
            return
    
    # Try to get the name
    try:
        metadata = identifiers.search_chemical(cas)
        if metadata and metadata.common_name:
            name = metadata.common_name
        else:
            name = cas_or_name
    except:
        name = cas_or_name
            
    print(f"\nAll Properties for: {name} (CAS: {cas})")
    print("=" * 60)
    
    # Dictionary of all property calculations
    properties = {
        "Critical Temperature": 
            (lambda: critical.Tc(cas), "K"),
        "Critical Pressure": 
            (lambda: critical.Pc(cas), "Pa"),
        "Critical Volume": 
            (lambda: critical.Vc(cas), "m³/mol"),
        "Critical Compressibility": 
            (lambda: critical.Zc(cas), "dimensionless"),
        "Boiling Point": 
            (lambda: phase_change.Tb(cas), "K"),
        "Melting Point": 
            (lambda: phase_change.Tm(cas), "K"),
        "Heat of Fusion": 
            (lambda: phase_change.Hfus(cas), "J/mol"),
        "Surface Tension": 
            (lambda: interface.sigma(cas, T=298.15), "N/m"),
        "Vapor Pressure": 
            (lambda: vapor_pressure.Psat(cas, T=298.15), "Pa"),
        "Acentric Factor": 
            (lambda: acentric.omega(cas), "dimensionless"),
        "Dipole Moment": 
            (lambda: dipole.dipole_moment(cas), "debye"),
        "Molecular Weight": 
            (lambda: identifiers.MW(cas), "g/mol"),
        "Heat of Formation (Gas)": 
            (lambda: reaction.Hfg(cas), "J/mol"),
        "Heat of Formation (Liquid)": 
            (lambda: reaction.Hfl(cas), "J/mol"),
        "Standard Entropy (Gas)": 
            (lambda: reaction.S0g(cas), "J/(mol·K)"),
        "Liquid Volume (COSTALD)": 
            (lambda: volume.COSTALD(T=298.15, Tc=critical.Tc(cas), Vc=critical.Vc(cas), omega=acentric.omega(cas)), "m³/mol"),
        "Liquid Density (Rackett)": 
            (lambda: 1.0/(volume.Rackett(T=298.15, Tc=critical.Tc(cas), Pc=critical.Pc(cas), omega=acentric.omega(cas))*identifiers.MW(cas)), "kg/m³"),
        "Gas Viscosity (Lucas)": 
            (lambda: viscosity.Lucas_gas(T=298.15, Tc=critical.Tc(cas), Pc=critical.Pc(cas), Zc=critical.Zc(cas), MW=identifiers.MW(cas), dipole=dipole.dipole_moment(cas)), "Pa·s"),
        "Liquid Viscosity (Letsou-Stiel)": 
            (lambda: viscosity.Letsou_Stiel(T=298.15, MW=identifiers.MW(cas), Tc=critical.Tc(cas), Pc=critical.Pc(cas), omega=acentric.omega(cas)), "Pa·s"),
        "Liquid Thermal Conductivity": 
            (lambda: thermal_conductivity.Sheffy_Johnson(T=298.15, MW=identifiers.MW(cas), Tc=critical.Tc(cas), Pc=critical.Pc(cas), omega=acentric.omega(cas)), "W/(m·K)"),
        "Gas Thermal Conductivity": 
            (lambda: thermal_conductivity.Eli_Hanley(T=298.15, MW=identifiers.MW(cas), Tc=critical.Tc(cas), Vc=critical.Vc(cas), Zc=critical.Zc(cas), omega=acentric.omega(cas), dipole=dipole.dipole_moment(cas)), "W/(m·K)"),
        "Hansen Solubility Parameter (Dispersion)": 
            (lambda: solubility.hansen_delta_d(cas), "Pa^0.5"),
        "Hansen Solubility Parameter (Polar)": 
            (lambda: solubility.hansen_delta_p(cas), "Pa^0.5"),
        "Hansen Solubility Parameter (Hydrogen Bonding)": 
            (lambda: solubility.hansen_delta_h(cas), "Pa^0.5"),
        "Triple Point Temperature": 
            (lambda: triple.Tt(cas), "K"),
        "Triple Point Pressure": 
            (lambda: triple.Pt(cas), "Pa"),
        "Flash Point": 
            (lambda: safety.Tflash(cas), "K"),
        "Autoignition Temperature": 
            (lambda: safety.Tautoignition(cas), "K"),
        "Lower Flammability Limit": 
            (lambda: safety.LFL(cas), "fraction"),
        "Upper Flammability Limit": 
            (lambda: safety.UFL(cas), "fraction"),
        "Global Warming Potential": 
            (lambda: environment.GWP(cas), "relative to CO₂"),
        "Ozone Depletion Potential": 
            (lambda: environment.ODP(cas), "relative to CFC-11"),
        "Log P (Octanol-Water Partition)": 
            (lambda: environment.logP(cas), "logarithmic")
    }
    
    # Calculate polarity classification separately
    try:
        dipole_value = dipole.dipole_moment(cas)
        if dipole_value is not None:
            polarity = determine_polarity(cas)
            properties["Polarity Classification"] = (lambda: polarity, "")
    except:
        pass
    
    # Calculate each property and display results
    results = []
    
    # First pass: get all properties
    for prop_name, (calc_func, unit) in properties.items():
        try:
            value = calc_func()
            if value is not None:
                results.append((prop_name, value, unit))
        except:
            pass
    
    # Sort results and display
    results.sort(key=lambda x: x[0])  # Sort alphabetically by property name
    
    # Group by categories for better organization
    categories = {
        "Critical Properties": ["Critical Temperature", "Critical Pressure", "Critical Volume", "Critical Compressibility"],
        "Phase Change": ["Boiling Point", "Melting Point", "Triple Point Temperature", "Triple Point Pressure"],
        "Thermodynamic": ["Heat of Fusion", "Heat of Formation (Gas)", "Heat of Formation (Liquid)", "Standard Entropy (Gas)"],
        "Transport Properties": ["Surface Tension", "Vapor Pressure", "Gas Viscosity", "Liquid Viscosity", "Gas Thermal Conductivity", "Liquid Thermal Conductivity"],
        "Physical Properties": ["Molecular Weight", "Acentric Factor", "Dipole Moment", "Polarity Classification", "Liquid Volume", "Liquid Density", "Hansen Solubility Parameter"],
        "Safety & Environmental": ["Flash Point", "Autoignition Temperature", "Lower Flammability Limit", "Upper Flammability Limit", "Global Warming Potential", "Ozone Depletion Potential", "Log P"]
    }
    
    # Display by category
    displayed_props = set()
    for category, props in categories.items():
        category_results = [r for r in results if any(p in r[0] for p in props)]
        if category_results:
            print(f"\n{category.upper()}:")
            for prop_name, value, unit in category_results:
                displayed_props.add(prop_name)
                if isinstance(value, float):
                    if abs(value) < 0.001 or abs(value) > 10000:
                        formatted_value = f"{value:.6e}"
                    else:
                        formatted_value = f"{value:.6g}"
                else:
                    formatted_value = str(value)
                
                if unit:
                    print(f"  {prop_name}: {formatted_value} {unit}")
                else:
                    print(f"  {prop_name}: {formatted_value}")
    
    # Display any remaining properties not in categories
    remaining = [r for r in results if r[0] not in displayed_props]
    if remaining:
        print("\nOTHER PROPERTIES:")
        for prop_name, value, unit in remaining:
            if isinstance(value, float):
                if abs(value) < 0.001 or abs(value) > 10000:
                    formatted_value = f"{value:.6e}"
                else:
                    formatted_value = f"{value:.6g}"
            else:
                formatted_value = str(value)
            
            if unit:
                print(f"  {prop_name}: {formatted_value} {unit}")
            else:
                print(f"  {prop_name}: {formatted_value}")
    
    # Show count of available properties
    print(f"\nTotal properties found: {len(results)} out of {len(properties)}")