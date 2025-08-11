#!/usr/bin/env python3
"""
Menu system for the Chemical Data CLI.
"""
import sys
from chemicals import identifiers

from chemdata.utils.validators import is_valid_cas
from chemdata.core.chemical_data import lookup_chemical, search_chemicals, get_all_properties
from chemdata.calculators.property_calculator import get_property_calculators, get_comparison_calculators, calculate_property, calculate_vapor_pressure_table
from chemdata.calculators.reaction_calculator import get_reaction_calculators, calculate_reaction_property
from chemdata.calculators.element_calculator import get_element_data
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

def vapor_pressure_table_menu():
    """Menu for calculating vapor pressure tables"""
    print_section_header("Vapor Pressure Table Calculator", "-")
    
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
            print("=" * 80)
            print(f"{'Temp (K)':>10} {'Temp (°C)':>10} {'Pressure (Pa)':>15} {'Pressure (bar)':>15} {'Pressure (mmHg)':>15} {'Pressure (kPa)':>15}")
            print("-" * 80)
            
            for row in results:
                print(f"{row['temperature_K']:>10.2f} {row['temperature_C']:>10.2f} "
                      f"{row['pressure_Pa']:>15.2e} {row['pressure_bar']:>15.6f} "
                      f"{row['pressure_mmHg']:>15.2f} {row['pressure_kPa']:>15.2f}")
            
            print(f"\nGenerated {len(results)} data points")
            
            # Ask if user wants to save to file
            save_choice = input("\nWould you like to save this table to a CSV file? (y/n): ").lower()
            if save_choice in ['y', 'yes']:
                filename = input("Enter filename (without extension): ") or f"vapor_pressure_{cas}"
                try:
                    import csv
                    with open(f"{filename}.csv", 'w', newline='') as csvfile:
                        fieldnames = ['temperature_K', 'temperature_C', 'pressure_Pa', 'pressure_bar', 'pressure_mmHg', 'pressure_kPa']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in results:
                            writer.writerow(row)
                    print(f"Table saved as {filename}.csv")
                except Exception as e:
                    print(f"Error saving file: {e}")
        else:
            print(f"No vapor pressure data available for {cas_or_name} in the specified temperature range.")
            print("This could be because:")
            print("- The chemical is not in the vapor pressure database")
            print("- The temperature range is outside the valid range for this chemical")
            print("- The chemical may be a solid at these temperatures")

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
9. Exit
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
            get_all_properties(cas_or_name)
            input("\nPress Enter to continue...")
            
        elif choice == '8':
            print("Exiting. Thank you for using Chemical Data CLI!")
            break
            
        else:
            print("Invalid choice, please try again.")