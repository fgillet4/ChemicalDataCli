#!/usr/bin/env python3
"""
Environmental impact analysis functionality for chemical compounds using chemicals.environment module.
"""
import sys

try:
    from chemicals import environment, identifiers
except ImportError:
    print("Required libraries not found. Make sure chemicals is installed.")
    print("Try: pip install chemicals thermo fluids ht")
    sys.exit(1)

from chemdata.utils.validators import is_valid_cas


def get_environmental_calculators():
    """Return a dictionary of environmental analysis functions with metadata."""
    return {
        "1": ("Global Warming Potential (GWP)", 
              lambda cas: environment.GWP(cas), 
              "[(impact/mass chemical)/(impact/mass CO2)]"),
        "2": ("Global Temperature Potential (GTP)", 
              lambda cas: environment.GTP(cas), 
              "[(impact/mass chemical)/(impact/mass CO2)]"),
        "3": ("Ozone Depletion Potential (ODP)", 
              lambda cas: environment.ODP(cas), 
              "[(impact/mass chemical)/(impact/mass CFC-11)]"),
        "4": ("Octanol-Water Partition Coefficient (logP)", 
              lambda cas: environment.logP(cas), 
              "log10(dimensionless)"),
        "5": ("GWP with specific method", 
              lambda cas: get_gwp_with_method(cas), 
              "[(impact/mass chemical)/(impact/mass CO2)]"),
        "6": ("GTP with specific method", 
              lambda cas: get_gtp_with_method(cas), 
              "[(impact/mass chemical)/(impact/mass CO2)]"),
        "7": ("ODP with specific method", 
              lambda cas: get_odp_with_method(cas), 
              "[(impact/mass chemical)/(impact/mass CFC-11)]"),
        "8": ("logP with specific method", 
              lambda cas: get_logp_with_method(cas), 
              "log10(dimensionless)")
    }


def get_gwp_with_method(cas):
    """Get GWP allowing user to select method."""
    available_methods = environment.GWP_methods(cas)
    if not available_methods:
        return "No methods available"
    
    print(f"\nAvailable GWP methods for {cas}:")
    for i, method in enumerate(available_methods, 1):
        print(f"{i}. {method}")
    
    try:
        choice = int(input("Select method (number): ")) - 1
        if 0 <= choice < len(available_methods):
            selected_method = available_methods[choice]
            return environment.GWP(cas, method=selected_method)
        else:
            print("Invalid selection, using default method")
            return environment.GWP(cas)
    except:
        print("Invalid input, using default method")
        return environment.GWP(cas)


def get_gtp_with_method(cas):
    """Get GTP allowing user to select method."""
    available_methods = environment.GTP_methods(cas)
    if not available_methods:
        return "No methods available"
    
    print(f"\nAvailable GTP methods for {cas}:")
    for i, method in enumerate(available_methods, 1):
        print(f"{i}. {method}")
    
    try:
        choice = int(input("Select method (number): ")) - 1
        if 0 <= choice < len(available_methods):
            selected_method = available_methods[choice]
            return environment.GTP(cas, method=selected_method)
        else:
            print("Invalid selection, using default method")
            return environment.GTP(cas)
    except:
        print("Invalid input, using default method")
        return environment.GTP(cas)


def get_odp_with_method(cas):
    """Get ODP allowing user to select method."""
    available_methods = environment.ODP_methods(cas)
    if not available_methods:
        return "No methods available"
    
    print(f"\nAvailable ODP methods for {cas}:")
    for i, method in enumerate(available_methods, 1):
        print(f"{i}. {method}")
    
    try:
        choice = int(input("Select method (number): ")) - 1
        if 0 <= choice < len(available_methods):
            selected_method = available_methods[choice]
            return environment.ODP(cas, method=selected_method)
        else:
            print("Invalid selection, using default method")
            return environment.ODP(cas)
    except:
        print("Invalid input, using default method")
        return environment.ODP(cas)


def get_logp_with_method(cas):
    """Get logP allowing user to select method."""
    available_methods = environment.logP_methods(cas)
    if not available_methods:
        return "No methods available"
    
    print(f"\nAvailable logP methods for {cas}:")
    for i, method in enumerate(available_methods, 1):
        print(f"{i}. {method}")
    
    try:
        choice = int(input("Select method (number): ")) - 1
        if 0 <= choice < len(available_methods):
            selected_method = available_methods[choice]
            return environment.logP(cas, method=selected_method)
        else:
            print("Invalid selection, using default method")
            return environment.logP(cas)
    except:
        print("Invalid input, using default method")
        return environment.logP(cas)


def calculate_environmental_property(property_key, cas_or_name):
    """Calculate a specific environmental property for a chemical."""
    env_props = get_environmental_calculators()
    
    if property_key not in env_props:
        print(f"Invalid environmental property key: {property_key}")
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
        property_name, calculation_func, unit_description = env_props[property_key]
        result = calculation_func(cas)
        
        return {
            "cas": cas,
            "name": cas_or_name,
            "property_name": property_name, 
            "value": result,
            "unit_description": unit_description
        }
    except Exception as e:
        print(f"Error calculating environmental property: {e}")
        return None


def comprehensive_environmental_analysis(cas_or_name):
    """Perform a comprehensive environmental impact analysis for a chemical."""
    # Validate and convert to CAS if necessary
    cas = None
    chemical_name = cas_or_name
    
    if is_valid_cas(cas_or_name):
        cas = cas_or_name
        try:
            chemical_name = identifiers.name(cas)
        except:
            chemical_name = cas_or_name
    else:
        try:
            cas = identifiers.CAS_from_any(cas_or_name)
        except Exception as e:
            print(f"Could not find chemical: {cas_or_name}, error: {e}")
            return None
    
    print(f"\n{'='*60}")
    print(f"COMPREHENSIVE ENVIRONMENTAL IMPACT ANALYSIS")
    print(f"{'='*60}")
    print(f"Chemical: {chemical_name}")
    print(f"CAS Number: {cas}")
    print(f"{'='*60}")
    
    # Initialize results
    results = {
        "cas": cas,
        "name": chemical_name,
        "climate_impact": {},
        "ozone_impact": {},
        "bioaccumulation": {},
        "environmental_recommendations": []
    }
    
    # 1. CLIMATE IMPACT
    print(f"\n1. CLIMATE IMPACT")
    print(f"{'-'*20}")
    
    # Global Warming Potential
    try:
        gwp_methods = environment.GWP_methods(cas)
        if gwp_methods:
            results["climate_impact"]["gwp_methods_available"] = gwp_methods
            
            # Get default GWP
            gwp_default = environment.GWP(cas)
            if gwp_default is not None:
                results["climate_impact"]["gwp_default"] = gwp_default
                print(f"Global Warming Potential (default): {gwp_default}")
                
                # Interpret GWP value
                if gwp_default > 1000:
                    print("  ⚠️  EXTREMELY HIGH climate impact")
                    results["environmental_recommendations"].append("CRITICAL: Extremely high GWP - avoid releases to atmosphere")
                elif gwp_default > 100:
                    print("  ⚠️  HIGH climate impact")
                    results["environmental_recommendations"].append("HIGH PRIORITY: High GWP - minimize atmospheric emissions")
                elif gwp_default > 10:
                    print("  ⚠️  MODERATE climate impact")
                    results["environmental_recommendations"].append("MODERATE: Consider climate impact in usage and disposal")
                elif gwp_default > 1:
                    print("  ⚠️  LOW climate impact")
                else:
                    print("  ✅ Minimal climate impact")
                
                # Get multiple timeframes if available
                for method in gwp_methods[:3]:  # Show first 3 methods
                    try:
                        gwp_specific = environment.GWP(cas, method=method)
                        if gwp_specific is not None:
                            results["climate_impact"][f"gwp_{method}"] = gwp_specific
                            print(f"  {method}: {gwp_specific}")
                    except:
                        pass
            else:
                print("Global Warming Potential: Not available")
        else:
            print("Global Warming Potential: No methods available")
    except Exception as e:
        print(f"Global Warming Potential: Error - {e}")
    
    # Global Temperature Potential
    try:
        gtp_methods = environment.GTP_methods(cas)
        if gtp_methods:
            results["climate_impact"]["gtp_methods_available"] = gtp_methods
            
            gtp_default = environment.GTP(cas)
            if gtp_default is not None:
                results["climate_impact"]["gtp_default"] = gtp_default
                print(f"Global Temperature Potential (default): {gtp_default}")
                
                # Show different timeframes
                for method in gtp_methods[:3]:
                    try:
                        gtp_specific = environment.GTP(cas, method=method)
                        if gtp_specific is not None:
                            results["climate_impact"][f"gtp_{method}"] = gtp_specific
                            print(f"  {method}: {gtp_specific}")
                    except:
                        pass
            else:
                print("Global Temperature Potential: Not available")
        else:
            print("Global Temperature Potential: No methods available")
    except Exception as e:
        print(f"Global Temperature Potential: Error - {e}")
    
    # 2. OZONE LAYER IMPACT
    print(f"\n2. OZONE LAYER IMPACT")
    print(f"{'-'*20}")
    
    try:
        odp_methods = environment.ODP_methods(cas)
        if odp_methods:
            results["ozone_impact"]["odp_methods_available"] = odp_methods
            
            odp_default = environment.ODP(cas)
            if odp_default is not None:
                results["ozone_impact"]["odp_default"] = odp_default
                print(f"Ozone Depletion Potential (default): {odp_default}")
                
                # Interpret ODP value
                if isinstance(odp_default, (int, float)):
                    if odp_default > 1.0:
                        print("  ⚠️  EXTREMELY HIGH ozone depletion risk")
                        results["environmental_recommendations"].append("CRITICAL: High ODP - substance may be regulated under Montreal Protocol")
                    elif odp_default > 0.1:
                        print("  ⚠️  HIGH ozone depletion risk")
                        results["environmental_recommendations"].append("HIGH PRIORITY: Significant ODP - avoid atmospheric releases")
                    elif odp_default > 0.01:
                        print("  ⚠️  MODERATE ozone depletion risk")
                        results["environmental_recommendations"].append("MODERATE: Monitor for ozone-depleting potential")
                    elif odp_default > 0:
                        print("  ⚠️  LOW ozone depletion risk")
                    else:
                        print("  ✅ No ozone depletion risk")
                
                # Show different methods
                for method in odp_methods[:3]:
                    try:
                        odp_specific = environment.ODP(cas, method=method)
                        if odp_specific is not None:
                            results["ozone_impact"][f"odp_{method}"] = odp_specific
                            print(f"  {method}: {odp_specific}")
                    except:
                        pass
            else:
                print("Ozone Depletion Potential: Not available")
        else:
            print("Ozone Depletion Potential: No methods available")
    except Exception as e:
        print(f"Ozone Depletion Potential: Error - {e}")
    
    # 3. BIOACCUMULATION POTENTIAL
    print(f"\n3. BIOACCUMULATION POTENTIAL")
    print(f"{'-'*20}")
    
    try:
        logp_methods = environment.logP_methods(cas)
        if logp_methods:
            results["bioaccumulation"]["logp_methods_available"] = logp_methods
            
            logp_default = environment.logP(cas)
            if logp_default is not None:
                results["bioaccumulation"]["logp_default"] = logp_default
                print(f"Octanol-Water Partition Coefficient (logP): {logp_default}")
                
                # Interpret logP value
                if logp_default > 5:
                    print("  ⚠️  VERY HIGH bioaccumulation potential")
                    results["environmental_recommendations"].append("CRITICAL: Very high bioaccumulation potential - environmental persistence concern")
                elif logp_default > 3:
                    print("  ⚠️  HIGH bioaccumulation potential")
                    results["environmental_recommendations"].append("HIGH PRIORITY: High bioaccumulation - monitor environmental releases")
                elif logp_default > 1:
                    print("  ⚠️  MODERATE bioaccumulation potential")
                    results["environmental_recommendations"].append("MODERATE: Moderate bioaccumulation potential")
                elif logp_default > -1:
                    print("  ✅ LOW bioaccumulation potential")
                else:
                    print("  ✅ VERY LOW bioaccumulation potential")
                
                # Additional interpretation
                if logp_default > 0:
                    print(f"  → More soluble in octanol (lipophilic)")
                else:
                    print(f"  → More soluble in water (hydrophilic)")
                
                # Show different data sources
                for method in logp_methods:
                    try:
                        logp_specific = environment.logP(cas, method=method)
                        if logp_specific is not None:
                            results["bioaccumulation"][f"logp_{method}"] = logp_specific
                            print(f"  {method}: {logp_specific}")
                    except:
                        pass
            else:
                print("Octanol-Water Partition Coefficient: Not available")
        else:
            print("Octanol-Water Partition Coefficient: No methods available")
    except Exception as e:
        print(f"Octanol-Water Partition Coefficient: Error - {e}")
    
    # 4. ENVIRONMENTAL RECOMMENDATIONS
    print(f"\n4. ENVIRONMENTAL RECOMMENDATIONS")
    print(f"{'-'*20}")
    
    if results["environmental_recommendations"]:
        for i, recommendation in enumerate(results["environmental_recommendations"], 1):
            print(f"{i}. {recommendation}")
    else:
        print("No specific environmental concerns identified from available data.")
    
    # General environmental recommendations
    print(f"\n5. GENERAL ENVIRONMENTAL PRACTICES")
    print(f"{'-'*20}")
    print("1. Follow all applicable environmental regulations")
    print("2. Implement proper waste management and disposal procedures")
    print("3. Monitor and minimize environmental releases")
    print("4. Consider green chemistry alternatives when possible")
    print("5. Maintain environmental incident response procedures")
    print("6. Regular environmental impact assessments")
    print("7. Employee training on environmental responsibilities")
    
    # 6. REGULATORY CONTEXT
    print(f"\n6. REGULATORY CONTEXT")
    print(f"{'-'*20}")
    
    # Check for potential regulatory concerns
    climate_concern = False
    ozone_concern = False
    
    if "gwp_default" in results["climate_impact"]:
        gwp = results["climate_impact"]["gwp_default"]
        if isinstance(gwp, (int, float)) and gwp > 100:
            climate_concern = True
    
    if "odp_default" in results["ozone_impact"]:
        odp = results["ozone_impact"]["odp_default"]
        if isinstance(odp, (int, float)) and odp > 0.01:
            ozone_concern = True
    
    if climate_concern:
        print("⚠️  May be subject to climate change regulations (high GWP)")
        print("  - Consider fluorinated gas regulations (F-gas)")
        print("  - May require emissions reporting")
    
    if ozone_concern:
        print("⚠️  May be subject to ozone protection regulations")
        print("  - Check Montreal Protocol restrictions")
        print("  - May require phase-out timeline compliance")
    
    if not climate_concern and not ozone_concern:
        print("✅ No obvious regulatory red flags identified")
        print("  - Always verify current local regulations")
        print("  - Environmental regulations change frequently")
    
    print(f"\n{'='*60}")
    print("END OF ENVIRONMENTAL ANALYSIS")
    print(f"{'='*60}")
    
    return results


def available_environmental_methods_for_chemical(cas_or_name):
    """Show which environmental analysis methods are available for a specific chemical."""
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
    
    print(f"\nAVAILABLE ENVIRONMENTAL METHODS FOR CAS: {cas}")
    print("=" * 50)
    
    methods = {
        "GWP (Global Warming Potential)": environment.GWP_methods,
        "GTP (Global Temperature Potential)": environment.GTP_methods,
        "ODP (Ozone Depletion Potential)": environment.ODP_methods,
        "logP (Octanol-Water Partition Coeff)": environment.logP_methods
    }
    
    for prop_name, method_func in methods.items():
        try:
            available = method_func(cas)
            
            if available:
                print(f"{prop_name}:")
                for method in available:
                    print(f"  - {method}")
            else:
                print(f"{prop_name}: No methods available")
            print()
        except Exception as e:
            print(f"{prop_name}: Error checking methods - {e}")
            print()


def environmental_comparison(chemicals_list):
    """Compare environmental impacts of multiple chemicals."""
    print(f"\nENVIRONMENTAL IMPACT COMPARISON")
    print("=" * 60)
    
    comparison_data = []
    
    for chem in chemicals_list:
        # Validate and convert to CAS if necessary
        cas = None
        if is_valid_cas(chem):
            cas = chem
        else:
            try:
                cas = identifiers.CAS_from_any(chem)
            except Exception as e:
                print(f"Could not find chemical: {chem}, error: {e}")
                continue
        
        # Get environmental properties
        chem_data = {"name": chem, "cas": cas}
        
        try:
            chem_data["gwp"] = environment.GWP(cas)
        except:
            chem_data["gwp"] = None
            
        try:
            chem_data["gtp"] = environment.GTP(cas)
        except:
            chem_data["gtp"] = None
            
        try:
            chem_data["odp"] = environment.ODP(cas)
        except:
            chem_data["odp"] = None
            
        try:
            chem_data["logp"] = environment.logP(cas)
        except:
            chem_data["logp"] = None
        
        comparison_data.append(chem_data)
    
    if not comparison_data:
        print("No valid chemicals to compare")
        return None
    
    # Display comparison table
    print(f"{'Chemical':>20} {'CAS':>15} {'GWP':>10} {'GTP':>10} {'ODP':>10} {'logP':>10}")
    print("-" * 80)
    
    for data in comparison_data:
        gwp_str = f"{data['gwp']:.1f}" if data['gwp'] is not None else "N/A"
        gtp_str = f"{data['gtp']:.1f}" if data['gtp'] is not None else "N/A"
        odp_str = f"{data['odp']:.3f}" if data['odp'] is not None and isinstance(data['odp'], (int, float)) else "N/A"
        logp_str = f"{data['logp']:.2f}" if data['logp'] is not None else "N/A"
        
        print(f"{data['name'][:20]:>20} {data['cas']:>15} {gwp_str:>10} {gtp_str:>10} {odp_str:>10} {logp_str:>10}")
    
    return comparison_data