#!/usr/bin/env python3
"""
Menu system for the Chemical Data CLI.
"""
import sys
from chemicals import identifiers

from chemdata.utils.validators import is_valid_cas
from chemdata.core.chemical_data import lookup_chemical, search_chemicals, get_all_properties
from chemdata.calculators.property_calculator import get_property_calculators, get_comparison_calculators, calculate_property, calculate_vapor_pressure_table, compare_vapor_pressure_methods, calculate_rotavapor_conditions, create_output_directory
from chemdata.calculators.reaction_calculator import get_reaction_calculators, calculate_reaction_property
from chemdata.calculators.element_calculator import get_element_data
from chemdata.calculators.safety_calculator import get_safety_calculators, calculate_safety_property, comprehensive_safety_analysis, available_methods_for_chemical, unit_conversion_helper
from chemdata.calculators.environmental_calculator import get_environmental_calculators, calculate_environmental_property, comprehensive_environmental_analysis, available_environmental_methods_for_chemical, environmental_comparison
from chemdata.ui.formatters import print_section_header, print_property, format_value

def property_calculator_menu():
    """Menu for calculating specific properties for any chemical"""
    properties = get_property_calculators()
    
    while True:
        print_section_header("Property Calculator", "-")
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
                    print(f"\n{property_name} of {cas_or_name} (CAS: {cas}): {format_value(result)} {unit}")
                    
                    # Additional information for dipole moment
                    if choice == "8":
                        from chemdata.core.chemical_data import determine_polarity
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
    comparison_calcs = get_comparison_calculators()
    
    while True:
        print_section_header("Chemical Property Comparison", "-")
        for key, calculator in comparison_calcs.items():
            print(f"{key}. Compare {calculator['name']}")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nEnter your choice: ").lower()
        
        if choice == 'back':
            break
            
        if choice in comparison_calcs:
            # Get the property function based on choice
            calculator = comparison_calcs[choice]
            prop_func = calculator['function']
            prop_name = calculator['name']
                
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
                        print(f"{name} (CAS: {cas}): {format_value(value)}")
                    else:
                        print(f"{name} (CAS: {cas}): Data not available")
            else:
                print("No valid results to compare")
        else:
            print("Invalid choice, please try again.")

def reaction_calculator_menu():
    """Calculate properties for chemical reactions"""
    reaction_calcs = get_reaction_calculators()
    
    while True:
        print_section_header("Chemical Reaction Calculator", "-")
        for key, calculator in reaction_calcs.items():
            print(f"{key}. Calculate {calculator['name']}")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nEnter your choice: ").lower()
        
        if choice == 'back':
            break
            
        if choice in reaction_calcs:
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
            T = float(input("Enter temperature in K (default is 298.15): ") or 298.15)
            
            try:
                # Calculate the property based on choice
                result = calculate_reaction_property(choice, cas, T)
                if result:
                    print(f"\n{result['property_name']} for {cas_or_name} at {T} K: {format_value(result['value'])} {result['unit']}")
            except Exception as e:
                print(f"Error calculating property: {e}")
        else:
            print("Invalid choice, please try again.")

def element_lookup_menu():
    """Menu for looking up element data from the periodic table"""
    print_section_header("Periodic Table Element Lookup", "-")
    
    while True:
        element = input("\nEnter element symbol, name, or atomic number (or 'back' to return): ")
        if element.lower() == 'back':
            break
            
        try:
            # Get element info
            element_data = get_element_data(element)
                
            if element_data:
                print(f"\nElement Information: {element_data['name']} ({element_data['symbol']})")
                print("=" * 50)
                print(f"Atomic Number: {element_data['atomic_number']}")
                print(f"Atomic Weight: {format_value(element_data['atomic_weight'])} g/mol")
                print(f"Period: {element_data['period']}")
                print(f"Group: {element_data['group']}")
                
                # Print block information
                if "block" in element_data:
                    print(f"Block: {element_data['block']}")
                
                # Print other available properties
                if "phase_at_stp" in element_data:
                    print(f"Phase at STP: {element_data['phase_at_stp']}")
                    
                if "heat_of_formation" in element_data:
                    print(f"Heat of Formation: {format_value(element_data['heat_of_formation'])} J/mol")
                    
                if "standard_entropy" in element_data:
                    print(f"Standard Entropy: {format_value(element_data['standard_entropy'])} J/(mol·K)")
                    
                if "electronegativity" in element_data:
                    print(f"Electronegativity: {format_value(element_data['electronegativity'])}")
                    
                if "covalent_radius" in element_data:
                    print(f"Covalent Radius: {format_value(element_data['covalent_radius'])} Å")
                    
                if "vdw_radius" in element_data:
                    print(f"Van der Waals Radius: {format_value(element_data['vdw_radius'])} Å")
            else:
                print(f"Could not find element: {element}")
                
        except Exception as e:
            print(f"Error looking up element {element}: {e}")

def safety_analysis_menu():
    """Menu for chemical safety analysis"""
    while True:
        print_section_header("Chemical Safety Analysis", "-")
        print("1. Individual safety property lookup")
        print("2. Comprehensive safety analysis")
        print("3. Available safety analysis methods for a chemical")
        print("4. Unit conversion helper")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nSelect option: ").lower()
        
        if choice == 'back':
            break
        elif choice == '1':
            safety_property_submenu()
        elif choice == '2':
            comprehensive_safety_submenu()
        elif choice == '3':
            available_methods_submenu()
        elif choice == '4':
            unit_conversion_helper()
        else:
            print("Invalid choice, please try again.")

def safety_property_submenu():
    """Submenu for individual safety property calculations"""
    safety_props = get_safety_calculators()
    
    while True:
        print_section_header("Individual Safety Properties", "-")
        for key, (name, _, unit_desc) in safety_props.items():
            print(f"{key}. {name}")
        print("Type 'back' to return to safety menu.")
        
        choice = input("\nSelect safety property to calculate: ").lower()
        
        if choice == 'back':
            break
            
        if choice in safety_props:
            cas_or_name = input("Enter chemical name or CAS number: ")
            
            result = calculate_safety_property(choice, cas_or_name)
            
            if result:
                property_name = result['property_name']
                value = result['value']
                unit_desc = result['unit_description']
                
                print(f"\n{property_name} for {result['name']} (CAS: {result['cas']}):")
                print("=" * 60)
                
                if isinstance(value, tuple) and len(value) == 2:
                    # For STEL, TWA, Ceiling which return (value, unit)
                    print(f"Value: {value[0]} {value[1]}")
                elif isinstance(value, bool):
                    # For skin absorption
                    print(f"Value: {'Yes' if value else 'No'}")
                elif isinstance(value, dict):
                    # For carcinogen status
                    print("Carcinogen Classifications:")
                    for source, status in value.items():
                        print(f"  - {source}: {status}")
                elif value is not None:
                    print(f"Value: {value}")
                    if unit_desc != "Returns (value, unit) where unit is ppm or mg/m³":
                        print(f"Unit: {unit_desc}")
                else:
                    print("No data available")
                
                print(f"Unit Description: {unit_desc}")
        else:
            print("Invalid choice, please try again.")

def comprehensive_safety_submenu():
    """Submenu for comprehensive safety analysis"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
        print("\nPerforming comprehensive safety analysis...")
        result = comprehensive_safety_analysis(cas_or_name)
        
        if result:
            # Ask if user wants to save to file
            save_choice = input("\nWould you like to save this safety analysis to a file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or f"safety_analysis_{result['cas']}"
                try:
                    import json
                    import os
                    from datetime import datetime
                    
                    # Create organized directory structure
                    cas_dir = create_output_directory(result['cas'])
                    
                    # Add timestamp to filename
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    full_filename = f"{filename}_{timestamp}.json"
                    filepath = os.path.join(cas_dir, full_filename)
                    
                    with open(filepath, 'w') as jsonfile:
                        json.dump(result, jsonfile, indent=2)
                    print(f"Safety analysis saved as: {filepath}")
                except Exception as e:
                    print(f"Error saving file: {e}")

def available_methods_submenu():
    """Submenu for showing available safety analysis methods"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
        available_methods_for_chemical(cas_or_name)

def environmental_analysis_menu():
    """Menu for environmental impact analysis"""
    while True:
        print_section_header("Environmental Impact Analysis", "-")
        print("1. Individual environmental property lookup")
        print("2. Comprehensive environmental analysis")
        print("3. Available environmental analysis methods for a chemical")
        print("4. Compare environmental impacts of multiple chemicals")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nSelect option: ").lower()
        
        if choice == 'back':
            break
        elif choice == '1':
            environmental_property_submenu()
        elif choice == '2':
            comprehensive_environmental_submenu()
        elif choice == '3':
            available_environmental_methods_submenu()
        elif choice == '4':
            environmental_comparison_submenu()
        else:
            print("Invalid choice, please try again.")

def environmental_property_submenu():
    """Submenu for individual environmental property calculations"""
    env_props = get_environmental_calculators()
    
    while True:
        print_section_header("Individual Environmental Properties", "-")
        for key, (name, _, unit_desc) in env_props.items():
            print(f"{key}. {name}")
        print("Type 'back' to return to environmental menu.")
        
        choice = input("\nSelect environmental property to calculate: ").lower()
        
        if choice == 'back':
            break
            
        if choice in env_props:
            cas_or_name = input("Enter chemical name or CAS number: ")
            
            result = calculate_environmental_property(choice, cas_or_name)
            
            if result:
                property_name = result['property_name']
                value = result['value']
                unit_desc = result['unit_description']
                
                print(f"\n{property_name} for {result['name']} (CAS: {result['cas']}):")
                print("=" * 60)
                
                if value is not None:
                    print(f"Value: {value}")
                    print(f"Unit: {unit_desc}")
                    
                    # Additional interpretation for specific properties
                    if "GWP" in property_name and isinstance(value, (int, float)):
                        if value > 1000:
                            print("Impact: EXTREMELY HIGH climate impact")
                        elif value > 100:
                            print("Impact: HIGH climate impact")
                        elif value > 10:
                            print("Impact: MODERATE climate impact")
                        elif value > 1:
                            print("Impact: LOW climate impact")
                        else:
                            print("Impact: Minimal climate impact")
                    
                    elif "ODP" in property_name and isinstance(value, (int, float)):
                        if value > 1.0:
                            print("Impact: EXTREMELY HIGH ozone depletion risk")
                        elif value > 0.1:
                            print("Impact: HIGH ozone depletion risk")
                        elif value > 0.01:
                            print("Impact: MODERATE ozone depletion risk")
                        elif value > 0:
                            print("Impact: LOW ozone depletion risk")
                        else:
                            print("Impact: No ozone depletion risk")
                    
                    elif "logP" in property_name and isinstance(value, (int, float)):
                        if value > 5:
                            print("Bioaccumulation: VERY HIGH potential")
                        elif value > 3:
                            print("Bioaccumulation: HIGH potential")
                        elif value > 1:
                            print("Bioaccumulation: MODERATE potential")
                        else:
                            print("Bioaccumulation: LOW potential")
                        
                        if value > 0:
                            print("Solubility: More lipophilic (fat-soluble)")
                        else:
                            print("Solubility: More hydrophilic (water-soluble)")
                else:
                    print("No data available")
        else:
            print("Invalid choice, please try again.")

def comprehensive_environmental_submenu():
    """Submenu for comprehensive environmental analysis"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
        print("\nPerforming comprehensive environmental analysis...")
        result = comprehensive_environmental_analysis(cas_or_name)
        
        if result:
            # Ask if user wants to save to file
            save_choice = input("\nWould you like to save this environmental analysis to a file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or f"environmental_analysis_{result['cas']}"
                try:
                    import json
                    import os
                    from datetime import datetime
                    
                    # Create organized directory structure
                    cas_dir = create_output_directory(result['cas'])
                    
                    # Add timestamp to filename
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    full_filename = f"{filename}_{timestamp}.json"
                    filepath = os.path.join(cas_dir, full_filename)
                    
                    with open(filepath, 'w') as jsonfile:
                        json.dump(result, jsonfile, indent=2)
                    print(f"Environmental analysis saved as: {filepath}")
                except Exception as e:
                    print(f"Error saving file: {e}")

def available_environmental_methods_submenu():
    """Submenu for showing available environmental analysis methods"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
        available_environmental_methods_for_chemical(cas_or_name)

def environmental_comparison_submenu():
    """Submenu for comparing environmental impacts of multiple chemicals"""
    while True:
        chemical_input = input("\nEnter chemical names or CAS numbers separated by commas (or 'back' to return): ")
        if chemical_input.lower() == 'back':
            break
            
        chemicals = [c.strip() for c in chemical_input.split(',')]
        
        if len(chemicals) < 2:
            print("Please enter at least 2 chemicals for comparison.")
            continue
        
        print("\nPerforming environmental comparison...")
        result = environmental_comparison(chemicals)
        
        if result:
            # Ask if user wants to save to file
            save_choice = input("\nWould you like to save this comparison to a CSV file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or "environmental_comparison"
                try:
                    import csv
                    import os
                    from datetime import datetime
                    
                    # Create directory if needed
                    if not os.path.exists("generated_data"):
                        os.makedirs("generated_data")
                    
                    # Add timestamp to filename
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    full_filename = f"{filename}_{timestamp}.csv"
                    filepath = os.path.join("generated_data", full_filename)
                    
                    with open(filepath, 'w', newline='') as csvfile:
                        fieldnames = ['name', 'cas', 'gwp', 'gtp', 'odp', 'logp']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in result:
                            writer.writerow(row)
                    print(f"Comparison saved as: {filepath}")
                except Exception as e:
                    print(f"Error saving file: {e}")

def vapor_pressure_menu():
    """Enhanced vapor pressure calculation menu with multiple options"""
    while True:
        print_section_header("Vapor Pressure Calculator", "-")
        print("1. Generate vapor pressure table over temperature range")
        print("2. Compare all vapor pressure calculation methods at single temperature")
        print("3. Rotavapor conditions (temperature/pressure combinations)")
        print("Type 'back' to return to the main menu.")
        
        choice = input("\nSelect option: ").lower()
        
        if choice == 'back':
            break
        elif choice == '1':
            vapor_pressure_table_submenu()
        elif choice == '2':
            vapor_pressure_comparison_submenu()
        elif choice == '3':
            rotavapor_submenu()
        else:
            print("Invalid choice, please try again.")

def vapor_pressure_table_submenu():
    """Submenu for generating vapor pressure tables"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
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
        
        # Get temperature range inputs
        try:
            T_start = float(input("Enter starting temperature in K (default 273.15): ") or 273.15)
            T_end = float(input("Enter ending temperature in K (default 373.15): ") or 373.15)
            T_step = float(input("Enter temperature step in K (default 10): ") or 10)
            
            if T_start >= T_end:
                print("Starting temperature must be less than ending temperature.")
                continue
                
            if T_step <= 0:
                print("Temperature step must be positive.")
                continue
                
        except ValueError:
            print("Invalid temperature input. Please enter numeric values.")
            continue
        
        # Calculate vapor pressure table
        print(f"\nCalculating vapor pressure table for {cas_or_name} (CAS: {cas})...")
        results = calculate_vapor_pressure_table(cas, T_start, T_end, T_step)
        
        if results and len(results) > 0:
            print(f"\nVapor Pressure Table for {cas_or_name}")
            print("=" * 100)
            print(f"{'Temp (K)':>10} {'Temp (°C)':>10} {'Pressure (Pa)':>15} {'Pressure (bar)':>15} {'Pressure (mmHg)':>15} {'Pressure (kPa)':>15} {'Method':>20}")
            print("-" * 100)
            
            for row in results:
                print(f"{row['temperature_K']:>10.2f} {row['temperature_C']:>10.2f} "
                      f"{row['pressure_Pa']:>15.2e} {row['pressure_bar']:>15.6f} "
                      f"{row['pressure_mmHg']:>15.2f} {row['pressure_kPa']:>15.2f} "
                      f"{row['method']:>20}")
            
            print(f"\nGenerated {len(results)} data points")
            
            # Ask if user wants to save to file
            save_choice = input("\nWould you like to save this table to a CSV file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or f"vapor_pressure_table_{T_start}K-{T_end}K"
                try:
                    import csv
                    import os
                    from datetime import datetime
                    
                    # Create organized directory structure
                    cas_dir = create_output_directory(cas)
                    
                    # Add timestamp to filename to avoid overwrites
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    full_filename = f"{filename}_{timestamp}.csv"
                    filepath = os.path.join(cas_dir, full_filename)
                    
                    with open(filepath, 'w', newline='') as csvfile:
                        fieldnames = ['temperature_K', 'temperature_C', 'pressure_Pa', 'pressure_bar', 'pressure_mmHg', 'pressure_kPa', 'method']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in results:
                            writer.writerow(row)
                    print(f"Table saved as: {filepath}")
                except Exception as e:
                    print(f"Error saving file: {e}")
        else:
            print(f"No vapor pressure data available for {cas_or_name} in the specified temperature range.")
            print("This could be because:")
            print("- The chemical is not in the vapor pressure database")
            print("- The temperature range is outside the valid range for this chemical")
            print("- The chemical may be a solid at these temperatures")

def vapor_pressure_comparison_submenu():
    """Submenu for comparing all vapor pressure calculation methods"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
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
        
        # Get temperature input
        try:
            T = float(input("Enter temperature in K (default 298.15): ") or 298.15)
        except ValueError:
            print("Invalid temperature input. Please enter a numeric value.")
            continue
        
        # Compare methods
        print(f"\nComparing vapor pressure methods for {cas_or_name} (CAS: {cas}) at {T} K...")
        results = compare_vapor_pressure_methods(cas, T)
        
        if results:
            print(f"\nVapor Pressure Method Comparison at {T} K ({T-273.15:.2f} °C)")
            print("=" * 80)
            print(f"{'Method':>25} {'Result (Pa)':>20} {'Result (mbar)':>15} {'Result (mmHg)':>15}")
            print("-" * 80)
            
            # Separate working methods from failed ones
            working_methods = {}
            failed_methods = {}
            
            for method, result in results.items():
                if isinstance(result, (int, float)) and result > 0:
                    working_methods[method] = result
                else:
                    failed_methods[method] = result
            
            # Display working methods first
            for method, pressure in working_methods.items():
                mbar = pressure / 100
                mmHg = pressure * 760 / 101325
                print(f"{method:>25} {pressure:>20.2e} {mbar:>15.2f} {mmHg:>15.2f}")
            
            # Display failed methods
            if failed_methods:
                print("\nFailed Methods:")
                print("-" * 80)
                for method, error in failed_methods.items():
                    print(f"{method:>25} {str(error):>50}")
            
            # Show statistics for working methods
            if len(working_methods) > 1:
                values = list(working_methods.values())
                avg = sum(values) / len(values)
                min_val = min(values)
                max_val = max(values)
                std_dev = (sum((x - avg) ** 2 for x in values) / len(values)) ** 0.5
                
                print(f"\nStatistics for working methods:")
                print(f"Average: {avg:.2e} Pa")
                print(f"Min: {min_val:.2e} Pa")
                print(f"Max: {max_val:.2e} Pa")
                print(f"Standard deviation: {std_dev:.2e} Pa")
                print(f"Coefficient of variation: {(std_dev/avg)*100:.1f}%")
            
            # Ask if user wants to save comparison to file
            save_choice = input("\nWould you like to save this comparison to a CSV file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or f"method_comparison_{T}K"
                try:
                    import csv
                    import os
                    from datetime import datetime
                    
                    # Create organized directory structure
                    cas_dir = create_output_directory(cas)
                    
                    # Add timestamp to filename
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    full_filename = f"{filename}_{timestamp}.csv"
                    filepath = os.path.join(cas_dir, full_filename)
                    
                    with open(filepath, 'w', newline='') as csvfile:
                        fieldnames = ['method', 'result_Pa', 'result_mbar', 'result_mmHg', 'status']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        
                        # Write working methods
                        for method, pressure in working_methods.items():
                            writer.writerow({
                                'method': method,
                                'result_Pa': pressure,
                                'result_mbar': pressure / 100,
                                'result_mmHg': pressure * 760 / 101325,
                                'status': 'Success'
                            })
                        
                        # Write failed methods
                        for method, error in failed_methods.items():
                            writer.writerow({
                                'method': method,
                                'result_Pa': 'N/A',
                                'result_mbar': 'N/A',
                                'result_mmHg': 'N/A',
                                'status': str(error)
                            })
                    
                    print(f"Comparison saved as: {filepath}")
                except Exception as e:
                    print(f"Error saving file: {e}")

def rotavapor_submenu():
    """Submenu for rotavapor conditions analysis"""
    while True:
        cas_or_name = input("\nEnter chemical name or CAS number (or 'back' to return): ")
        if cas_or_name.lower() == 'back':
            break
            
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
        
        print("\nRotavapor Condition Analysis")
        print("Enter temperature and pressure ranges for your rotavapor setup")
        
        try:
            # Temperature inputs
            T_input = input("Enter temperatures in K (e.g., '293,313,333' or '293-333-20' for range): ")
            if '-' in T_input:
                parts = T_input.split('-')
                T_start, T_end, T_step = float(parts[0]), float(parts[1]), float(parts[2])
                T_values = [T_start + i * T_step for i in range(int((T_end - T_start) / T_step) + 1)]
            else:
                T_values = [float(T.strip()) for T in T_input.split(',')]
            
            # Pressure inputs
            P_input = input("Enter pressures in Pa (e.g., '1000,5000,10000' or '1000-50000-5000' for range): ")
            if '-' in P_input:
                parts = P_input.split('-')
                P_start, P_end, P_step = float(parts[0]), float(parts[1]), float(parts[2])
                P_values = [P_start + i * P_step for i in range(int((P_end - P_start) / P_step) + 1)]
            else:
                P_values = [float(P.strip()) for P in P_input.split(',')]
                
        except ValueError:
            print("Invalid input format. Please use commas or range format (start-end-step).")
            continue
        
        # Calculate rotavapor conditions
        print(f"\nAnalyzing rotavapor conditions for {cas_or_name} (CAS: {cas})...")
        results = calculate_rotavapor_conditions(cas, T_values, P_values)
        
        if results and len(results) > 0:
            print(f"\nRotavapor Conditions Analysis for {cas_or_name}")
            print("=" * 130)
            print(f"{'T(K)':>6} {'T(°C)':>7} {'Sys P(Pa)':>10} {'Sys P(mbar)':>12} {'Vap P(Pa)':>12} {'Vap P(mbar)':>12} {'P Ratio':>8} {'Evaporate':>10} {'Rate':>10} {'Method':>15}")
            print("-" * 130)
            
            for row in results:
                evap_symbol = "✓" if row['will_evaporate'] else "✗"
                print(f"{row['temperature_K']:>6.1f} {row['temperature_C']:>7.1f} "
                      f"{row['system_pressure_Pa']:>10.0f} {row['system_pressure_mbar']:>12.1f} "
                      f"{row['vapor_pressure_Pa']:>12.2e} {row['vapor_pressure_mbar']:>12.1f} "
                      f"{row['pressure_ratio']:>8.2f} {evap_symbol:>10} "
                      f"{row['evaporation_rate']:>10} {row['method']:>15}")
            
            print(f"\nGenerated {len(results)} condition combinations")
            print("\nLegend:")
            print("✓ = Will evaporate (vapor pressure > system pressure)")
            print("✗ = Will not evaporate")
            print("Rate: None < Low < Moderate < High")
            
            # Ask if user wants to save to file
            save_choice = input("\nWould you like to save this analysis to a CSV file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or f"rotavapor_analysis"
                try:
                    import csv
                    import os
                    from datetime import datetime
                    
                    # Create organized directory structure
                    cas_dir = create_output_directory(cas)
                    
                    # Add timestamp and temperature/pressure info to filename
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    T_range = f"{min(T_values):.0f}K-{max(T_values):.0f}K"
                    P_range = f"{min(P_values):.0f}Pa-{max(P_values):.0f}Pa"
                    full_filename = f"{filename}_{T_range}_{P_range}_{timestamp}.csv"
                    filepath = os.path.join(cas_dir, full_filename)
                    
                    with open(filepath, 'w', newline='') as csvfile:
                        fieldnames = ['temperature_K', 'temperature_C', 'system_pressure_Pa', 'system_pressure_mbar',
                                     'vapor_pressure_Pa', 'vapor_pressure_mbar', 'pressure_ratio', 
                                     'will_evaporate', 'evaporation_rate', 'method']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in results:
                            writer.writerow(row)
                    print(f"Analysis saved as: {filepath}")
                except Exception as e:
                    print(f"Error saving file: {e}")
        else:
            print(f"No vapor pressure data available for {cas_or_name}.")
            print("Try using a different chemical or check the CAS number.")

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
7. Get ALL available properties for a chemical
8. Calculate vapor pressure tables
9. Chemical safety analysis
10. Environmental impact analysis
11. Exit
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
            reaction_calculator_menu()
            
        elif choice == '6':
            element_lookup_menu()
            
        elif choice == '7':
            cas_or_name = input("Enter chemical name or CAS number: ")
            temp_input = input("Enter temperature in K (press Enter for 298.15 K / 25°C): ").strip()
            if temp_input:
                try:
                    T = float(temp_input)
                    get_all_properties(cas_or_name, T=T)
                except ValueError:
                    print("Invalid temperature. Using default 298.15 K.")
                    get_all_properties(cas_or_name)
            else:
                get_all_properties(cas_or_name)
            input("\nPress Enter to continue...")
            
        elif choice == '8':
            vapor_pressure_menu()
            
        elif choice == '9':
            safety_analysis_menu()
            
        elif choice == '10':
            environmental_analysis_menu()
            
        elif choice == '11':
            print("Exiting. Thank you for using Chemical Data CLI!")
            break
            
        else:
            print("Invalid choice, please try again.")