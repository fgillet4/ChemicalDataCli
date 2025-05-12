#!/usr/bin/env python3
"""
Formatters for displaying chemical data in the CLI.
"""

def format_value(value):
    """Format a value for display, handling different data types appropriately."""
    if value is None:
        return "N/A"
    
    if isinstance(value, float):
        if abs(value) < 0.001 or abs(value) > 10000:
            return f"{value:.6e}"
        else:
            return f"{value:.6g}"
    elif isinstance(value, list) and len(value) > 5:
        return f"{value[:5]} ... ({len(value)} items)"
    else:
        return str(value)
        
def format_property_name(name):
    """Format a property name for display, converting snake_case to Title Case."""
    return " ".join(word.capitalize() for word in name.split("_"))

def print_section_header(title, char="=", width=50):
    """Print a section header with a title and underline."""
    print(f"\n{title}")
    print(char * width)
    
def print_property(name, value, unit=None):
    """Print a property with optional unit."""
    formatted_value = format_value(value)
    if unit:
        print(f"  {name}: {formatted_value} {unit}")
    else:
        print(f"  {name}: {formatted_value}")
        
def print_property_dict(properties, title=None):
    """Print a dictionary of properties with a title."""
    if title:
        print_section_header(title)
        
    for name, value in properties.items():
        if isinstance(value, dict) and "value" in value and "unit" in value:
            print_property(format_property_name(name), value["value"], value["unit"])
        else:
            print_property(format_property_name(name), value)
            
def print_categorized_properties(properties, categories, title=None):
    """Print properties organized by categories."""
    if title:
        print_section_header(title)
        
    # Track properties that have been displayed
    displayed_props = set()
    
    # Display properties by category
    for category_name, category_props in categories.items():
        # Get properties that belong to this category
        matching_props = {}
        for prop_name, prop_value in properties.items():
            if any(marker in prop_name for marker in category_props):
                matching_props[prop_name] = prop_value
                displayed_props.add(prop_name)
                
        if matching_props:
            print(f"\n{category_name.upper()}:")
            for name, value in matching_props.items():
                if isinstance(value, dict) and "value" in value and "unit" in value:
                    print_property(format_property_name(name), value["value"], value["unit"])
                else:
                    print_property(format_property_name(name), value)
    
    # Display any remaining properties
    remaining_props = {name: value for name, value in properties.items() 
                       if name not in displayed_props}
    if remaining_props:
        print("\nOTHER PROPERTIES:")
        for name, value in remaining_props.items():
            if isinstance(value, dict) and "value" in value and "unit" in value:
                print_property(format_property_name(name), value["value"], value["unit"])
            else:
                print_property(format_property_name(name), value)