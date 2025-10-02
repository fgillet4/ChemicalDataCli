#!/usr/bin/env python3
"""
Comprehensive thermodynamic properties calculator using the thermo library.
Includes heat capacity, thermal conductivity, viscosity, and other core properties.
"""
import sys

try:
    from thermo import heat_capacity, vapor_pressure, eos
    from chemicals import identifiers, critical
    import numpy as np
except ImportError:
    print("Required libraries not found. Make sure thermo and chemicals are installed.")
    print("Try: pip install thermo chemicals")
    sys.exit(1)

from chemdata.utils.validators import is_valid_cas


def get_thermo_property_calculators():
    """Return a dictionary of thermodynamic property calculators with metadata."""
    return {
        "1": ("Liquid Heat Capacity vs Temperature", 
              lambda cas, T: calculate_heat_capacity_liquid(cas, T), 
              "J/mol/K"),
        "2": ("Gas Heat Capacity vs Temperature", 
              lambda cas, T: calculate_heat_capacity_gas(cas, T), 
              "J/mol/K"),
        "3": ("Solid Heat Capacity vs Temperature", 
              lambda cas, T: calculate_heat_capacity_solid(cas, T), 
              "J/mol/K"),
        "4": ("Heat Capacity Temperature Profile", 
              lambda cas: generate_heat_capacity_profile(cas), 
              "Temperature profile data"),
        "5": ("Compare Heat Capacity Methods", 
              lambda cas, T, phase: compare_heat_capacity_methods(cas, T, phase), 
              "Method comparison"),
        "6": ("Heat Capacity Integration (Enthalpy Change)", 
              lambda cas, T1, T2, phase: integrate_heat_capacity(cas, T1, T2, phase), 
              "J/mol"),
        "7": ("Available Heat Capacity Methods", 
              lambda cas: get_available_heat_capacity_methods(cas), 
              "Method list"),
        "8": ("Vapor Pressure vs Temperature", 
              lambda cas, T: calculate_vapor_pressure(cas, T), 
              "Pa"),
        "9": ("Sublimation Pressure vs Temperature", 
              lambda cas, T: calculate_sublimation_pressure(cas, T), 
              "Pa"),
        "10": ("Vapor Pressure Temperature Profile", 
               lambda cas: generate_vapor_pressure_profile(cas), 
               "Temperature profile data"),
        "11": ("Compare Vapor Pressure Methods", 
               lambda cas, T: compare_vapor_pressure_methods(cas, T), 
               "Method comparison"),
        "12": ("Antoine Equation Parameters", 
               lambda cas: get_antoine_parameters(cas), 
               "Antoine coefficients"),
        "13": ("Boiling Point from Vapor Pressure", 
               lambda cas, P: calculate_boiling_point(cas, P), 
               "K"),
        "14": ("Available Vapor Pressure Methods", 
               lambda cas: get_available_vapor_pressure_methods(cas), 
               "Method list"),
        "15": ("Peng-Robinson EOS Properties", 
               lambda cas, T, P: calculate_eos_properties(cas, T, P, 'PR'), 
               "EOS state properties"),
        "16": ("Soave-Redlich-Kwong EOS Properties", 
               lambda cas, T, P: calculate_eos_properties(cas, T, P, 'SRK'), 
               "EOS state properties"),
        "17": ("Van der Waals EOS Properties", 
               lambda cas, T, P: calculate_eos_properties(cas, T, P, 'VDW'), 
               "EOS state properties"),
        "18": ("Compare Multiple EOS Models", 
               lambda cas, T, P: compare_eos_models(cas, T, P), 
               "EOS comparison"),
        "19": ("Saturation Properties from EOS", 
               lambda cas, T: calculate_saturation_properties(cas, T), 
               "Saturation state properties"),
        "20": ("Phase Identification", 
               lambda cas, T, P: identify_phase_eos(cas, T, P), 
               "Phase identification"),
        "21": ("Critical Point Calculations", 
               lambda cas: calculate_critical_properties(cas), 
               "Critical properties"),
        "22": ("EOS Volume Solutions", 
               lambda cas, T, P: analyze_eos_volumes(cas, T, P), 
               "Volume analysis")
    }


def calculate_heat_capacity_liquid(cas, T):
    """Calculate liquid heat capacity at specified temperature."""
    try:
        # Get molecular weight for the calculation
        MW = identifiers.MW(cas)
        
        # Create HeatCapacityLiquid object
        Cp_liquid = heat_capacity.HeatCapacityLiquid(CASRN=cas, MW=MW)
        
        # Calculate at specified temperature
        result = Cp_liquid(T)
        
        return {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "heat_capacity": result,
            "method_used": Cp_liquid.method,
            "available_methods": Cp_liquid.all_methods,
            "valid_temperature_range": {
                "Tmin": Cp_liquid.Tmin,
                "Tmax": Cp_liquid.Tmax
            }
        }
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_heat_capacity_gas(cas, T):
    """Calculate gas heat capacity at specified temperature."""
    try:
        # Get molecular weight for the calculation
        MW = identifiers.MW(cas)
        
        # Create HeatCapacityGas object
        Cp_gas = heat_capacity.HeatCapacityGas(CASRN=cas, MW=MW)
        
        # Calculate at specified temperature
        result = Cp_gas(T)
        
        return {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "heat_capacity": result,
            "method_used": Cp_gas.method,
            "available_methods": Cp_gas.all_methods,
            "valid_temperature_range": {
                "Tmin": Cp_gas.Tmin,
                "Tmax": Cp_gas.Tmax
            }
        }
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_heat_capacity_solid(cas, T):
    """Calculate solid heat capacity at specified temperature."""
    try:
        # Get molecular weight for the calculation
        MW = identifiers.MW(cas)
        
        # Create HeatCapacitySolid object
        Cp_solid = heat_capacity.HeatCapacitySolid(CASRN=cas, MW=MW)
        
        # Calculate at specified temperature
        result = Cp_solid(T)
        
        return {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "heat_capacity": result,
            "method_used": Cp_solid.method,
            "available_methods": Cp_solid.all_methods,
            "valid_temperature_range": {
                "Tmin": Cp_solid.Tmin,
                "Tmax": Cp_solid.Tmax
            }
        }
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def generate_heat_capacity_profile(cas, T_start=200, T_end=800, num_points=50):
    """Generate heat capacity vs temperature profiles for all phases."""
    try:
        # Get molecular weight
        MW = identifiers.MW(cas)
        
        # Create objects for all phases
        Cp_liquid = heat_capacity.HeatCapacityLiquid(CASRN=cas, MW=MW)
        Cp_gas = heat_capacity.HeatCapacityGas(CASRN=cas, MW=MW)
        Cp_solid = heat_capacity.HeatCapacitySolid(CASRN=cas, MW=MW)
        
        # Temperature range
        temperatures = np.linspace(T_start, T_end, num_points)
        
        results = {
            "temperatures_K": temperatures.tolist(),
            "temperatures_C": (temperatures - 273.15).tolist(),
            "liquid": {"values": [], "valid_range": None, "method": None},
            "gas": {"values": [], "valid_range": None, "method": None},
            "solid": {"values": [], "valid_range": None, "method": None}
        }
        
        # Calculate liquid heat capacities
        for T in temperatures:
            try:
                if Cp_liquid.Tmin <= T <= Cp_liquid.Tmax:
                    cp_val = Cp_liquid(T)
                    results["liquid"]["values"].append(cp_val)
                else:
                    results["liquid"]["values"].append(None)
            except:
                results["liquid"]["values"].append(None)
        
        results["liquid"]["valid_range"] = {"Tmin": Cp_liquid.Tmin, "Tmax": Cp_liquid.Tmax}
        results["liquid"]["method"] = Cp_liquid.method
        
        # Calculate gas heat capacities
        for T in temperatures:
            try:
                if Cp_gas.Tmin <= T <= Cp_gas.Tmax:
                    cp_val = Cp_gas(T)
                    results["gas"]["values"].append(cp_val)
                else:
                    results["gas"]["values"].append(None)
            except:
                results["gas"]["values"].append(None)
        
        results["gas"]["valid_range"] = {"Tmin": Cp_gas.Tmin, "Tmax": Cp_gas.Tmax}
        results["gas"]["method"] = Cp_gas.method
        
        # Calculate solid heat capacities
        for T in temperatures:
            try:
                if Cp_solid.Tmin <= T <= Cp_solid.Tmax:
                    cp_val = Cp_solid(T)
                    results["solid"]["values"].append(cp_val)
                else:
                    results["solid"]["values"].append(None)
            except:
                results["solid"]["values"].append(None)
        
        results["solid"]["valid_range"] = {"Tmin": Cp_solid.Tmin, "Tmax": Cp_solid.Tmax}
        results["solid"]["method"] = Cp_solid.method
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_heat_capacity_methods(cas, T, phase):
    """Compare all available heat capacity calculation methods at a given temperature."""
    try:
        MW = identifiers.MW(cas)
        
        if phase.lower() == 'liquid':
            cp_obj = heat_capacity.HeatCapacityLiquid(CASRN=cas, MW=MW)
        elif phase.lower() == 'gas':
            cp_obj = heat_capacity.HeatCapacityGas(CASRN=cas, MW=MW)
        elif phase.lower() == 'solid':
            cp_obj = heat_capacity.HeatCapacitySolid(CASRN=cas, MW=MW)
        else:
            return {"error": "Phase must be 'liquid', 'gas', or 'solid'"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "phase": phase,
            "methods": {},
            "default_method": cp_obj.method,
            "default_value": None
        }
        
        # Get default value
        try:
            results["default_value"] = cp_obj(T)
        except:
            results["default_value"] = None
        
        # Test all available methods
        for method in cp_obj.all_methods:
            try:
                # Test if method is valid at this temperature
                if cp_obj.test_method_validity(T, method):
                    value = cp_obj.calculate(T, method)
                    results["methods"][method] = {
                        "value": value,
                        "status": "success"
                    }
                else:
                    results["methods"][method] = {
                        "value": None,
                        "status": "invalid_temperature"
                    }
            except Exception as e:
                results["methods"][method] = {
                    "value": None,
                    "status": f"error: {str(e)}"
                }
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def integrate_heat_capacity(cas, T1, T2, phase, num_points=100):
    """Integrate heat capacity to calculate enthalpy change between two temperatures."""
    try:
        MW = identifiers.MW(cas)
        
        if phase.lower() == 'liquid':
            cp_obj = heat_capacity.HeatCapacityLiquid(CASRN=cas, MW=MW)
        elif phase.lower() == 'gas':
            cp_obj = heat_capacity.HeatCapacityGas(CASRN=cas, MW=MW)
        elif phase.lower() == 'solid':
            cp_obj = heat_capacity.HeatCapacitySolid(CASRN=cas, MW=MW)
        else:
            return {"error": "Phase must be 'liquid', 'gas', or 'solid'"}
        
        # Check temperature range validity
        T_min = max(T1, cp_obj.Tmin) if cp_obj.Tmin else T1
        T_max = min(T2, cp_obj.Tmax) if cp_obj.Tmax else T2
        
        if T_min >= T_max:
            return {"error": f"Invalid temperature range for this method. Valid range: {cp_obj.Tmin}-{cp_obj.Tmax} K"}
        
        # Numerical integration using trapezoidal rule
        temperatures = np.linspace(T_min, T_max, num_points)
        cp_values = []
        
        for T in temperatures:
            try:
                cp_val = cp_obj(T)
                cp_values.append(cp_val)
            except:
                return {"error": f"Failed to calculate heat capacity at {T} K"}
        
        cp_values = np.array(cp_values)
        
        # Integrate using trapezoidal rule
        enthalpy_change = np.trapz(cp_values, temperatures)
        
        results = {
            "T1_K": T1,
            "T2_K": T2,
            "T1_C": T1 - 273.15,
            "T2_C": T2 - 273.15,
            "phase": phase,
            "enthalpy_change_J_mol": enthalpy_change,
            "enthalpy_change_kJ_mol": enthalpy_change / 1000,
            "method_used": cp_obj.method,
            "integration_points": num_points,
            "valid_integration_range": {
                "T_min_K": T_min,
                "T_max_K": T_max
            }
        }
        
        # Add per-gram values if molecular weight is available
        if MW:
            results["enthalpy_change_J_g"] = enthalpy_change / MW
            results["enthalpy_change_kJ_kg"] = enthalpy_change / MW
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def get_available_heat_capacity_methods(cas):
    """Get all available heat capacity calculation methods for a chemical."""
    try:
        MW = identifiers.MW(cas)
        
        # Create objects for all phases
        Cp_liquid = heat_capacity.HeatCapacityLiquid(CASRN=cas, MW=MW)
        Cp_gas = heat_capacity.HeatCapacityGas(CASRN=cas, MW=MW)
        Cp_solid = heat_capacity.HeatCapacitySolid(CASRN=cas, MW=MW)
        
        results = {
            "liquid": {
                "available_methods": Cp_liquid.all_methods,
                "default_method": Cp_liquid.method,
                "temperature_range": {
                    "Tmin": Cp_liquid.Tmin,
                    "Tmax": Cp_liquid.Tmax
                }
            },
            "gas": {
                "available_methods": Cp_gas.all_methods,
                "default_method": Cp_gas.method,
                "temperature_range": {
                    "Tmin": Cp_gas.Tmin,
                    "Tmax": Cp_gas.Tmax
                }
            },
            "solid": {
                "available_methods": Cp_solid.all_methods,
                "default_method": Cp_solid.method,
                "temperature_range": {
                    "Tmin": Cp_solid.Tmin,
                    "Tmax": Cp_solid.Tmax
                }
            }
        }
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_thermo_property(property_key, cas_or_name, **kwargs):
    """Calculate a specific thermodynamic property for a chemical."""
    thermo_props = get_thermo_property_calculators()
    
    if property_key not in thermo_props:
        print(f"Invalid thermodynamic property key: {property_key}")
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
        property_name, calculation_func, unit_description = thermo_props[property_key]
        
        # Handle different argument requirements for different properties
        if property_key in ["1", "2", "3"]:  # Heat capacity at specific temperature
            T = kwargs.get('temperature', 298.15)
            result = calculation_func(cas, T)
        elif property_key == "4":  # Temperature profile
            result = calculation_func(cas)
        elif property_key == "5":  # Method comparison
            T = kwargs.get('temperature', 298.15)
            phase = kwargs.get('phase', 'gas')
            result = calculation_func(cas, T, phase)
        elif property_key == "6":  # Integration
            T1 = kwargs.get('T1', 298.15)
            T2 = kwargs.get('T2', 373.15)
            phase = kwargs.get('phase', 'gas')
            result = calculation_func(cas, T1, T2, phase)
        elif property_key == "7":  # Available methods
            result = calculation_func(cas)
        elif property_key in ["8", "9"]:  # Vapor/sublimation pressure at temperature
            T = kwargs.get('temperature', 298.15)
            result = calculation_func(cas, T)
        elif property_key == "10":  # Vapor pressure profile
            result = calculation_func(cas)
        elif property_key == "11":  # Method comparison
            T = kwargs.get('temperature', 298.15)
            result = calculation_func(cas, T)
        elif property_key in ["12", "14"]:  # Parameters or available methods
            result = calculation_func(cas)
        elif property_key == "13":  # Boiling point calculation
            P = kwargs.get('pressure', 101325)
            result = calculation_func(cas, P)
        elif property_key in ["15", "16", "17", "18", "20", "22"]:  # EOS calculations requiring T,P
            T = kwargs.get('temperature', 298.15)
            P = kwargs.get('pressure', 101325)
            result = calculation_func(cas, T, P)
        elif property_key == "19":  # Saturation properties
            T = kwargs.get('temperature', 298.15)
            result = calculation_func(cas, T)
        elif property_key == "21":  # Critical properties
            result = calculation_func(cas)
        else:
            result = calculation_func(cas)
        
        return {
            "cas": cas,
            "name": cas_or_name,
            "property_name": property_name,
            "result": result,
            "unit_description": unit_description
        }
    except Exception as e:
        print(f"Error calculating thermodynamic property: {e}")
        return None


# ====================== VAPOR PRESSURE CALCULATIONS ======================

def calculate_vapor_pressure(cas, T):
    """Calculate vapor pressure at specified temperature."""
    try:
        # Get critical properties for correlations
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        # Create VaporPressure object
        vp_obj = vapor_pressure.VaporPressure(CASRN=cas, Tc=Tc, Pc=Pc, omega=omega)
        
        # Calculate at specified temperature
        result = vp_obj(T)
        
        return {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "vapor_pressure_Pa": result,
            "vapor_pressure_bar": result / 100000 if result else None,
            "vapor_pressure_mmHg": result * 760 / 101325 if result else None,
            "vapor_pressure_kPa": result / 1000 if result else None,
            "method_used": vp_obj.method,
            "available_methods": vp_obj.all_methods,
            "valid_temperature_range": {
                "Tmin": vp_obj.Tmin,
                "Tmax": vp_obj.Tmax
            }
        }
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_sublimation_pressure(cas, T):
    """Calculate sublimation pressure at specified temperature."""
    try:
        # Get triple point properties if available
        from chemicals import triple
        try:
            Tt = triple.Tt(cas)
            Pt = triple.Pt(cas)
        except:
            Tt = None
            Pt = None
        
        # Create SublimationPressure object
        sub_obj = vapor_pressure.SublimationPressure(CASRN=cas, Tt=Tt, Pt=Pt)
        
        # Calculate at specified temperature
        result = sub_obj(T)
        
        return {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "sublimation_pressure_Pa": result,
            "sublimation_pressure_bar": result / 100000 if result else None,
            "sublimation_pressure_mmHg": result * 760 / 101325 if result else None,
            "sublimation_pressure_kPa": result / 1000 if result else None,
            "method_used": sub_obj.method,
            "available_methods": sub_obj.all_methods,
            "valid_temperature_range": {
                "Tmin": sub_obj.Tmin,
                "Tmax": sub_obj.Tmax
            }
        }
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def generate_vapor_pressure_profile(cas, T_start=250, T_end=600, num_points=50):
    """Generate vapor pressure vs temperature profile."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        # Create VaporPressure object
        vp_obj = vapor_pressure.VaporPressure(CASRN=cas, Tc=Tc, Pc=Pc, omega=omega)
        
        # Adjust temperature range based on valid range
        if vp_obj.Tmin and T_start < vp_obj.Tmin:
            T_start = vp_obj.Tmin
        if vp_obj.Tmax and T_end > vp_obj.Tmax:
            T_end = vp_obj.Tmax
        
        if T_start >= T_end:
            return {"error": f"Invalid temperature range. Valid range: {vp_obj.Tmin}-{vp_obj.Tmax} K"}
        
        # Temperature range
        temperatures = np.linspace(T_start, T_end, num_points)
        
        results = {
            "temperatures_K": temperatures.tolist(),
            "temperatures_C": (temperatures - 273.15).tolist(),
            "vapor_pressures_Pa": [],
            "vapor_pressures_bar": [],
            "vapor_pressures_mmHg": [],
            "vapor_pressures_kPa": [],
            "method_used": vp_obj.method,
            "valid_range": {
                "Tmin": vp_obj.Tmin,
                "Tmax": vp_obj.Tmax
            }
        }
        
        # Calculate vapor pressures
        for T in temperatures:
            try:
                vp = vp_obj(T)
                results["vapor_pressures_Pa"].append(vp)
                results["vapor_pressures_bar"].append(vp / 100000 if vp else None)
                results["vapor_pressures_mmHg"].append(vp * 760 / 101325 if vp else None)
                results["vapor_pressures_kPa"].append(vp / 1000 if vp else None)
            except:
                results["vapor_pressures_Pa"].append(None)
                results["vapor_pressures_bar"].append(None)
                results["vapor_pressures_mmHg"].append(None)
                results["vapor_pressures_kPa"].append(None)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_vapor_pressure_methods(cas, T):
    """Compare all available vapor pressure calculation methods at a given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        vp_obj = vapor_pressure.VaporPressure(CASRN=cas, Tc=Tc, Pc=Pc, omega=omega)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "methods": {},
            "default_method": vp_obj.method,
            "default_value": None
        }
        
        # Get default value
        try:
            results["default_value"] = vp_obj(T)
        except:
            results["default_value"] = None
        
        # Test all available methods
        for method in vp_obj.all_methods:
            try:
                # Test if method is valid at this temperature
                if vp_obj.test_method_validity(T, method):
                    value = vp_obj.calculate(T, method)
                    results["methods"][method] = {
                        "value_Pa": value,
                        "value_bar": value / 100000 if value else None,
                        "value_mmHg": value * 760 / 101325 if value else None,
                        "value_kPa": value / 1000 if value else None,
                        "status": "success"
                    }
                else:
                    results["methods"][method] = {
                        "value_Pa": None,
                        "value_bar": None,
                        "value_mmHg": None,
                        "value_kPa": None,
                        "status": "invalid_temperature"
                    }
            except Exception as e:
                results["methods"][method] = {
                    "value_Pa": None,
                    "value_bar": None,
                    "value_mmHg": None,
                    "value_kPa": None,
                    "status": f"error: {str(e)}"
                }
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def get_antoine_parameters(cas):
    """Get Antoine equation parameters for the chemical."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        vp_obj = vapor_pressure.VaporPressure(CASRN=cas, Tc=Tc, Pc=Pc, omega=omega)
        
        results = {
            "antoine_methods": [],
            "parameters": {}
        }
        
        # Check for Antoine-based methods
        antoine_methods = [method for method in vp_obj.all_methods if 'ANTOINE' in method]
        results["antoine_methods"] = antoine_methods
        
        # Try to get coefficients for each Antoine method
        for method in antoine_methods:
            try:
                # The specific coefficients would need to be extracted from the object's data
                # This is a simplified version - actual implementation would access internal data
                results["parameters"][method] = {
                    "note": "Antoine parameters available but extraction requires detailed implementation",
                    "form": "ln(P) = A - B/(T + C) where P is in Pa and T in K"
                }
            except Exception as e:
                results["parameters"][method] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_boiling_point(cas, P=101325):
    """Calculate boiling point at specified pressure using vapor pressure correlation."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        vp_obj = vapor_pressure.VaporPressure(CASRN=cas, Tc=Tc, Pc=Pc, omega=omega)
        
        # Use solve_property to find temperature at given pressure
        try:
            T_boiling = vp_obj.solve_property(P)
            
            return {
                "pressure_Pa": P,
                "pressure_bar": P / 100000,
                "pressure_mmHg": P * 760 / 101325,
                "pressure_kPa": P / 1000,
                "boiling_point_K": T_boiling,
                "boiling_point_C": T_boiling - 273.15,
                "method_used": vp_obj.method
            }
        except:
            # If solve_property fails, use iterative method
            T_guess = 373.15  # Start with 100°C
            tolerance = 1.0  # Pa
            max_iterations = 100
            
            for i in range(max_iterations):
                try:
                    P_calc = vp_obj(T_guess)
                    if abs(P_calc - P) < tolerance:
                        return {
                            "pressure_Pa": P,
                            "pressure_bar": P / 100000,
                            "pressure_mmHg": P * 760 / 101325,
                            "pressure_kPa": P / 1000,
                            "boiling_point_K": T_guess,
                            "boiling_point_C": T_guess - 273.15,
                            "method_used": vp_obj.method,
                            "iterations": i + 1
                        }
                    
                    # Simple Newton-like iteration
                    if P_calc > P:
                        T_guess -= 1.0
                    else:
                        T_guess += 1.0
                        
                except:
                    break
            
            return {"error": "Could not converge to boiling point"}
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def get_available_vapor_pressure_methods(cas):
    """Get all available vapor pressure calculation methods for a chemical."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        vp_obj = vapor_pressure.VaporPressure(CASRN=cas, Tc=Tc, Pc=Pc, omega=omega)
        sub_obj = vapor_pressure.SublimationPressure(CASRN=cas)
        
        results = {
            "vapor_pressure": {
                "available_methods": vp_obj.all_methods,
                "default_method": vp_obj.method,
                "temperature_range": {
                    "Tmin": vp_obj.Tmin,
                    "Tmax": vp_obj.Tmax
                }
            },
            "sublimation_pressure": {
                "available_methods": sub_obj.all_methods,
                "default_method": sub_obj.method,
                "temperature_range": {
                    "Tmin": sub_obj.Tmin,
                    "Tmax": sub_obj.Tmax
                }
            }
        }
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


# ====================== EOS CALCULATIONS ======================

def calculate_eos_properties(cas, T, P, eos_type='PR'):
    """Calculate thermodynamic properties using cubic equations of state."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Create EOS object based on type
        if eos_type.upper() == 'PR':
            eos_obj = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'SRK':
            eos_obj = eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        elif eos_type.upper() == 'VDW':
            eos_obj = eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        else:
            return {"error": f"Unsupported EOS type: {eos_type}"}
        
        # Calculate properties
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "eos_type": eos_type,
            "critical_properties": {
                "Tc_K": Tc,
                "Pc_Pa": Pc,
                "omega": omega
            }
        }
        
        # Volume solutions
        try:
            volumes = eos_obj.volume_solutions(T, P)
            if volumes:
                results["volume_solutions"] = {
                    "all_volumes_m3_mol": volumes,
                    "number_of_solutions": len(volumes)
                }
                
                # Identify phases for each volume
                valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
                if valid_volumes:
                    min_vol = min(valid_volumes)
                    max_vol = max(valid_volumes)
                    results["liquid_volume_m3_mol"] = min_vol
                    results["gas_volume_m3_mol"] = max_vol
        except Exception as e:
            results["volume_solutions"] = {"error": str(e)}
        
        # Fugacity
        try:
            fugacity_liquid = eos_obj.fugacity_l(T, P)
            fugacity_gas = eos_obj.fugacity_g(T, P)
            results["fugacity"] = {
                "liquid_Pa": fugacity_liquid,
                "gas_Pa": fugacity_gas,
                "fugacity_coefficient_liquid": fugacity_liquid / P if fugacity_liquid else None,
                "fugacity_coefficient_gas": fugacity_gas / P if fugacity_gas else None
            }
        except Exception as e:
            results["fugacity"] = {"error": str(e)}
        
        # Compressibility factor
        try:
            Z_l = eos_obj.Z_l
            Z_g = eos_obj.Z_g
            results["compressibility_factor"] = {
                "liquid": Z_l,
                "gas": Z_g
            }
        except Exception as e:
            results["compressibility_factor"] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def compare_eos_models(cas, T, P):
    """Compare multiple EOS models at given conditions."""
    try:
        eos_models = ['PR', 'SRK', 'VDW']
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "models": {}
        }
        
        for model in eos_models:
            try:
                model_result = calculate_eos_properties(cas, T, P, model)
                if "error" not in model_result:
                    results["models"][model] = {
                        "liquid_volume": model_result.get("liquid_volume_m3_mol"),
                        "gas_volume": model_result.get("gas_volume_m3_mol"),
                        "Z_liquid": model_result.get("compressibility_factor", {}).get("liquid"),
                        "Z_gas": model_result.get("compressibility_factor", {}).get("gas"),
                        "fugacity_liquid": model_result.get("fugacity", {}).get("liquid_Pa"),
                        "fugacity_gas": model_result.get("fugacity", {}).get("gas_Pa")
                    }
                else:
                    results["models"][model] = {"error": model_result["error"]}
            except Exception as e:
                results["models"][model] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_saturation_properties(cas, T):
    """Calculate saturation properties using EOS at given temperature."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        if T >= Tc:
            return {"error": f"Temperature {T} K is above critical temperature {Tc} K"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "critical_temperature_K": Tc
        }
        
        try:
            # Calculate saturation pressure
            Psat = pr_eos.Psat(T)
            results["saturation_pressure_Pa"] = Psat
            results["saturation_pressure_bar"] = Psat / 100000
            results["saturation_pressure_mmHg"] = Psat * 760 / 101325
            
            # Calculate saturated volumes
            V_l_sat = pr_eos.V_l_sat(T)
            V_g_sat = pr_eos.V_g_sat(T)
            results["saturated_liquid_volume_m3_mol"] = V_l_sat
            results["saturated_gas_volume_m3_mol"] = V_g_sat
            
            # Calculate compressibility factors
            Z_l_sat = Psat * V_l_sat / (8.314 * T)
            Z_g_sat = Psat * V_g_sat / (8.314 * T)
            results["Z_liquid_saturated"] = Z_l_sat
            results["Z_gas_saturated"] = Z_g_sat
            
        except Exception as e:
            results["calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def identify_phase_eos(cas, T, P):
    """Identify phase state using EOS at given conditions."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        # Use Peng-Robinson EOS
        pr_eos = eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P)
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "critical_temperature_K": Tc,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000
        }
        
        # Determine phase based on critical conditions
        if T > Tc and P > Pc:
            phase = "Supercritical fluid"
        elif T > Tc:
            phase = "Supercritical gas"
        else:
            try:
                # Calculate saturation pressure at this temperature
                Psat = pr_eos.Psat(T)
                results["saturation_pressure_Pa"] = Psat
                results["saturation_pressure_bar"] = Psat / 100000
                
                if P > Psat:
                    phase = "Liquid"
                elif P < Psat:
                    phase = "Gas"
                else:
                    phase = "Two-phase (vapor-liquid equilibrium)"
                    
            except Exception as e:
                # If saturation pressure calculation fails, use simple heuristics
                if P > Pc * 0.5:  # Rough estimate
                    phase = "Likely liquid"
                else:
                    phase = "Likely gas"
                results["saturation_calculation_error"] = str(e)
        
        results["identified_phase"] = phase
        
        # Calculate volume to verify phase identification
        try:
            volumes = pr_eos.volume_solutions(T, P)
            valid_volumes = [v for v in volumes if v > 0 and not np.isnan(v)]
            
            if valid_volumes:
                results["volume_solutions_m3_mol"] = valid_volumes
                if len(valid_volumes) == 1:
                    results["phase_verification"] = "Single-phase confirmed"
                elif len(valid_volumes) == 3:
                    results["phase_verification"] = "Three volume solutions - two-phase region possible"
                    
        except Exception as e:
            results["volume_calculation_error"] = str(e)
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def calculate_critical_properties(cas):
    """Calculate and analyze critical properties."""
    try:
        # Get critical properties from database
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        Vc = critical.Vc(cas)
        omega = critical.omega(cas)
        
        results = {
            "cas": cas,
            "critical_temperature_K": Tc,
            "critical_temperature_C": Tc - 273.15 if Tc else None,
            "critical_pressure_Pa": Pc,
            "critical_pressure_bar": Pc / 100000 if Pc else None,
            "critical_volume_m3_mol": Vc,
            "acentric_factor": omega
        }
        
        # Calculate derived properties
        if all([Tc, Pc, Vc]):
            # Critical compressibility factor
            Zc = Pc * Vc / (8.314 * Tc)
            results["critical_compressibility_factor"] = Zc
            
            # Reduced properties at standard conditions (298.15 K, 101325 Pa)
            T_std = 298.15
            P_std = 101325
            results["reduced_temperature_at_STP"] = T_std / Tc
            results["reduced_pressure_at_STP"] = P_std / Pc
            
        # Get available calculation methods
        try:
            methods = {
                "Tc_methods": critical.Tc_methods(cas),
                "Pc_methods": critical.Pc_methods(cas),
                "Vc_methods": critical.Vc_methods(cas),
                "omega_methods": critical.omega_methods(cas)
            }
            results["available_methods"] = methods
        except:
            pass
        
        return results
        
    except Exception as e:
        return {"error": str(e)}


def analyze_eos_volumes(cas, T, P):
    """Analyze volume solutions from cubic EOS."""
    try:
        # Get critical properties
        from chemicals import critical
        Tc = critical.Tc(cas)
        Pc = critical.Pc(cas)
        omega = critical.omega(cas)
        
        if not all([Tc, Pc, omega]):
            return {"error": "Critical properties not available for this chemical"}
        
        results = {
            "temperature_K": T,
            "temperature_C": T - 273.15,
            "pressure_Pa": P,
            "pressure_bar": P / 100000,
            "analysis": {}
        }
        
        # Analyze different EOS models
        eos_models = {
            'PR': eos.PR(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'SRK': eos.SRK(Tc=Tc, Pc=Pc, omega=omega, T=T, P=P),
            'VDW': eos.VDW(Tc=Tc, Pc=Pc, T=T, P=P)
        }
        
        for model_name, eos_obj in eos_models.items():
            try:
                # Get volume solutions
                volumes = eos_obj.volume_solutions(T, P)
                
                model_analysis = {
                    "raw_solutions": volumes,
                    "number_of_solutions": len(volumes),
                    "physical_solutions": [],
                    "complex_solutions": [],
                    "negative_solutions": []
                }
                
                # Classify solutions
                for vol in volumes:
                    if np.isnan(vol) or np.isinf(vol):
                        continue
                    elif np.iscomplex(vol):
                        model_analysis["complex_solutions"].append(complex(vol))
                    elif vol < 0:
                        model_analysis["negative_solutions"].append(vol)
                    else:
                        model_analysis["physical_solutions"].append(vol)
                
                # Identify liquid and gas phases from physical solutions
                physical_vols = model_analysis["physical_solutions"]
                if len(physical_vols) >= 2:
                    model_analysis["liquid_volume_m3_mol"] = min(physical_vols)
                    model_analysis["gas_volume_m3_mol"] = max(physical_vols)
                    model_analysis["phase_behavior"] = "Two-phase possible"
                elif len(physical_vols) == 1:
                    vol = physical_vols[0]
                    # Use critical volume to estimate phase
                    if hasattr(eos_obj, 'Vc') and eos_obj.Vc:
                        if vol < eos_obj.Vc:
                            model_analysis["phase_behavior"] = "Likely liquid"
                        else:
                            model_analysis["phase_behavior"] = "Likely gas"
                    model_analysis["single_volume_m3_mol"] = vol
                else:
                    model_analysis["phase_behavior"] = "No physical solutions"
                
                # Calculate compressibility factors for physical solutions
                if physical_vols:
                    Z_values = [P * v / (8.314 * T) for v in physical_vols]
                    model_analysis["compressibility_factors"] = Z_values
                
                results["analysis"][model_name] = model_analysis
                
            except Exception as e:
                results["analysis"][model_name] = {"error": str(e)}
        
        return results
        
    except Exception as e:
        return {"error": str(e)}