#!/usr/bin/env python3
"""
Safety analysis functionality for chemical compounds using chemicals.safety module.
"""
import sys

try:
    from chemicals import safety, identifiers
except ImportError:
    print("Required libraries not found. Make sure chemicals is installed.")
    print("Try: pip install chemicals thermo fluids ht")
    sys.exit(1)

from chemdata.utils.validators import is_valid_cas


def get_safety_calculators():
    """Return a dictionary of safety analysis functions with metadata."""
    return {
        "1": ("Short-term Exposure Limit (STEL)", 
              lambda cas: safety.STEL(cas), 
              "Returns (value, unit) where unit is ppm or mg/m³"),
        "2": ("Time-Weighted Average Exposure Limit (TWA)", 
              lambda cas: safety.TWA(cas), 
              "Returns (value, unit) where unit is ppm or mg/m³"),
        "3": ("Ceiling Exposure Limit", 
              lambda cas: safety.Ceiling(cas), 
              "Returns (value, unit) where unit is ppm or mg/m³"),
        "4": ("Skin Absorption", 
              lambda cas: safety.Skin(cas), 
              "Returns True/False if absorbed through skin"),
        "5": ("Carcinogen Status", 
              lambda cas: safety.Carcinogen(cas), 
              "Returns carcinogen classification from multiple sources"),
        "6": ("Flash Point", 
              lambda cas: safety.T_flash(cas), 
              "K"),
        "7": ("Autoignition Temperature", 
              lambda cas: safety.T_autoignition(cas), 
              "K"),
        "8": ("Lower Flammability Limit (LFL)", 
              lambda cas: safety.LFL(CASRN=cas), 
              "mole fraction"),
        "9": ("Upper Flammability Limit (UFL)", 
              lambda cas: safety.UFL(CASRN=cas), 
              "mole fraction"),
        "10": ("NFPA 30 Classification", 
               lambda cas: classify_nfpa_30(cas), 
               "NFPA 30 flammability class")
    }


def classify_nfpa_30(cas):
    """Classify chemical according to NFPA 30 flammability standard."""
    try:
        # Get flash point and boiling point
        T_flash = safety.T_flash(cas)
        if T_flash is None:
            return "Data not available"
        
        # Try to get boiling point for better classification
        from chemicals import phase_change
        Tb = phase_change.Tb(cas)
        
        # Use NFPA 30 classification
        classification = safety.NFPA_30_classification(T_flash, Tb=Tb)
        return classification
        
    except Exception as e:
        return f"Classification failed: {str(e)}"


def calculate_safety_property(property_key, cas_or_name):
    """Calculate a specific safety property for a chemical."""
    safety_props = get_safety_calculators()
    
    if property_key not in safety_props:
        print(f"Invalid safety property key: {property_key}")
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
        property_name, calculation_func, unit_description = safety_props[property_key]
        result = calculation_func(cas)
        
        return {
            "cas": cas,
            "name": cas_or_name,
            "property_name": property_name, 
            "value": result,
            "unit_description": unit_description
        }
    except Exception as e:
        print(f"Error calculating safety property: {e}")
        return None


def comprehensive_safety_analysis(cas_or_name):
    """Perform a comprehensive safety analysis for a chemical."""
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
    print(f"COMPREHENSIVE SAFETY ANALYSIS")
    print(f"{'='*60}")
    print(f"Chemical: {chemical_name}")
    print(f"CAS Number: {cas}")
    print(f"{'='*60}")
    
    # Initialize results
    results = {
        "cas": cas,
        "name": chemical_name,
        "exposure_limits": {},
        "physical_hazards": {},
        "health_hazards": {},
        "flammability": {},
        "recommendations": []
    }
    
    # 1. EXPOSURE LIMITS
    print(f"\n1. EXPOSURE LIMITS")
    print(f"{'-'*20}")
    
    # STEL
    try:
        stel_result = safety.STEL(cas)
        if stel_result:
            value, unit = stel_result
            results["exposure_limits"]["STEL"] = {"value": value, "unit": unit}
            print(f"Short-term Exposure Limit (STEL): {value} {unit}")
        else:
            print("Short-term Exposure Limit (STEL): Not available")
    except Exception as e:
        print(f"Short-term Exposure Limit (STEL): Error - {e}")
    
    # TWA
    try:
        twa_result = safety.TWA(cas)
        if twa_result:
            value, unit = twa_result
            results["exposure_limits"]["TWA"] = {"value": value, "unit": unit}
            print(f"Time-Weighted Average (TWA): {value} {unit}")
        else:
            print("Time-Weighted Average (TWA): Not available")
    except Exception as e:
        print(f"Time-Weighted Average (TWA): Error - {e}")
    
    # Ceiling
    try:
        ceiling_result = safety.Ceiling(cas)
        if ceiling_result:
            value, unit = ceiling_result
            results["exposure_limits"]["Ceiling"] = {"value": value, "unit": unit}
            print(f"Ceiling Limit: {value} {unit}")
        else:
            print("Ceiling Limit: Not available")
    except Exception as e:
        print(f"Ceiling Limit: Error - {e}")
    
    # 2. HEALTH HAZARDS
    print(f"\n2. HEALTH HAZARDS")
    print(f"{'-'*20}")
    
    # Skin absorption
    try:
        skin_absorption = safety.Skin(cas)
        if skin_absorption is not None:
            results["health_hazards"]["skin_absorption"] = skin_absorption
            print(f"Skin Absorption: {'Yes' if skin_absorption else 'No'}")
            if skin_absorption:
                results["recommendations"].append("Use appropriate skin protection (gloves, protective clothing)")
        else:
            print("Skin Absorption: Data not available")
    except Exception as e:
        print(f"Skin Absorption: Error - {e}")
    
    # Carcinogenicity
    try:
        carcinogen_data = safety.Carcinogen(cas)
        if carcinogen_data:
            results["health_hazards"]["carcinogenicity"] = carcinogen_data
            print("Carcinogenicity Status:")
            if isinstance(carcinogen_data, dict):
                for source, status in carcinogen_data.items():
                    print(f"  - {source}: {status}")
                    if any(keyword in status.lower() for keyword in ['carcinogen', 'cancer', 'reasonably anticipated']):
                        results["recommendations"].append("Handle as potential carcinogen - use fume hood and minimize exposure")
            else:
                print(f"  - {carcinogen_data}")
        else:
            print("Carcinogenicity Status: Not available")
    except Exception as e:
        print(f"Carcinogenicity Status: Error - {e}")
    
    # 3. FLAMMABILITY HAZARDS
    print(f"\n3. FLAMMABILITY HAZARDS")
    print(f"{'-'*20}")
    
    # Flash Point
    try:
        flash_point = safety.T_flash(cas)
        if flash_point:
            flash_point_c = flash_point - 273.15
            results["flammability"]["flash_point"] = {"K": flash_point, "C": flash_point_c}
            print(f"Flash Point: {flash_point:.1f} K ({flash_point_c:.1f} °C)")
            
            # Flash point safety recommendations
            if flash_point_c < 23:
                results["recommendations"].append("EXTREMELY FLAMMABLE - Keep away from heat, sparks, and open flames")
            elif flash_point_c < 60:
                results["recommendations"].append("FLAMMABLE - Keep away from ignition sources")
            elif flash_point_c < 93:
                results["recommendations"].append("COMBUSTIBLE - Avoid heating near ignition sources")
        else:
            print("Flash Point: Not available")
    except Exception as e:
        print(f"Flash Point: Error - {e}")
    
    # Autoignition Temperature
    try:
        autoignition_temp = safety.T_autoignition(cas)
        if autoignition_temp:
            autoignition_c = autoignition_temp - 273.15
            results["flammability"]["autoignition_temp"] = {"K": autoignition_temp, "C": autoignition_c}
            print(f"Autoignition Temperature: {autoignition_temp:.1f} K ({autoignition_c:.1f} °C)")
        else:
            print("Autoignition Temperature: Not available")
    except Exception as e:
        print(f"Autoignition Temperature: Error - {e}")
    
    # Lower Flammability Limit
    try:
        lfl = safety.LFL(CASRN=cas)
        if lfl:
            lfl_percent = lfl * 100
            results["flammability"]["LFL"] = {"fraction": lfl, "percent": lfl_percent}
            print(f"Lower Flammability Limit: {lfl:.4f} ({lfl_percent:.2f}%)")
        else:
            print("Lower Flammability Limit: Not available")
    except Exception as e:
        print(f"Lower Flammability Limit: Error - {e}")
    
    # Upper Flammability Limit
    try:
        ufl = safety.UFL(CASRN=cas)
        if ufl:
            ufl_percent = ufl * 100
            results["flammability"]["UFL"] = {"fraction": ufl, "percent": ufl_percent}
            print(f"Upper Flammability Limit: {ufl:.4f} ({ufl_percent:.2f}%)")
        else:
            print("Upper Flammability Limit: Not available")
    except Exception as e:
        print(f"Upper Flammability Limit: Error - {e}")
    
    # NFPA 30 Classification
    try:
        nfpa_class = classify_nfpa_30(cas)
        if nfpa_class and "not available" not in nfpa_class.lower():
            results["flammability"]["NFPA_30_class"] = nfpa_class
            print(f"NFPA 30 Classification: Class {nfpa_class}")
            
            # NFPA class specific recommendations
            if nfpa_class in ['IA', 'IB', 'IC']:
                results["recommendations"].append(f"NFPA Class {nfpa_class} FLAMMABLE LIQUID - Store in approved flammable storage")
            elif nfpa_class in ['II', 'IIIA', 'IIIB']:
                results["recommendations"].append(f"NFPA Class {nfpa_class} COMBUSTIBLE LIQUID - Store away from ignition sources")
        else:
            print("NFPA 30 Classification: Not available")
    except Exception as e:
        print(f"NFPA 30 Classification: Error - {e}")
    
    # 4. SAFETY RECOMMENDATIONS
    print(f"\n4. SAFETY RECOMMENDATIONS")
    print(f"{'-'*20}")
    
    if results["recommendations"]:
        for i, recommendation in enumerate(results["recommendations"], 1):
            print(f"{i}. {recommendation}")
    else:
        print("No specific recommendations available based on current data.")
        print("Always follow general laboratory safety practices:")
        print("1. Use appropriate personal protective equipment")
        print("2. Work in well-ventilated areas")
        print("3. Keep containers tightly closed")
        print("4. Wash hands thoroughly after handling")
    
    # General recommendations
    print(f"\n5. GENERAL SAFETY PRACTICES")
    print(f"{'-'*20}")
    print("1. Read and understand the Safety Data Sheet (SDS)")
    print("2. Use appropriate personal protective equipment (PPE)")
    print("3. Ensure adequate ventilation or use fume hood")
    print("4. Keep emergency contact information readily available")
    print("5. Know the location of safety equipment (eyewash, shower, fire extinguisher)")
    print("6. Store according to compatibility guidelines")
    print("7. Dispose of waste according to institutional guidelines")
    
    print(f"\n{'='*60}")
    print("END OF SAFETY ANALYSIS")
    print(f"{'='*60}")
    
    return results


def unit_conversion_helper():
    """Provide unit conversion utilities for exposure limits."""
    print("\nUNIT CONVERSION HELPER")
    print("=" * 30)
    print("For exposure limit conversions:")
    print("- ppm to mg/m³: Use ppmv_to_mgm3(ppmv, MW, T, P)")
    print("- mg/m³ to ppm: Use mgm3_to_ppmv(mgm3, MW, T, P)")
    print("- Default conditions: T=298.15 K, P=101325 Pa")
    print("\nExample:")
    print("from chemicals.safety import ppmv_to_mgm3, mgm3_to_ppmv")
    print("ppmv_to_mgm3(1.0, 40.0)  # 1 ppm of MW=40 g/mol gas")
    print("mgm3_to_ppmv(1.635, 40.0)  # 1.635 mg/m³ of MW=40 g/mol gas")


def available_methods_for_chemical(cas_or_name):
    """Show which safety analysis methods are available for a specific chemical."""
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
    
    print(f"\nAVAILABLE SAFETY METHODS FOR CAS: {cas}")
    print("=" * 50)
    
    methods = {
        "STEL": safety.STEL_methods,
        "TWA": safety.TWA_methods,
        "Ceiling": safety.Ceiling_methods,
        "Skin": safety.Skin_methods,
        "Carcinogen": safety.Carcinogen_methods,
        "Flash Point": safety.T_flash_methods,
        "Autoignition": safety.T_autoignition_methods,
        "LFL": safety.LFL_methods,
        "UFL": safety.UFL_methods
    }
    
    for prop_name, method_func in methods.items():
        try:
            if prop_name in ["LFL", "UFL"]:
                available = method_func(CASRN=cas)
            else:
                available = method_func(cas)
            
            if available:
                print(f"{prop_name}: {', '.join(available)}")
            else:
                print(f"{prop_name}: No methods available")
        except Exception as e:
            print(f"{prop_name}: Error checking methods - {e}")