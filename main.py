#!/usr/bin/env python3
"""
Chemical Data CLI - A tool for retrieving and calculating chemical properties.

This application provides a command-line interface for accessing chemical property data 
from the chemicals library.
"""
import sys
import argparse

try:
    from chemdata.ui.menus import main_menu
    from chemdata.core.chemical_data import lookup_chemical, search_chemicals, get_all_properties
    from chemdata import __version__
except ImportError:
    print("Error importing required modules. Make sure the 'chemdata' package is in your Python path.")
    sys.exit(1)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Chemical Data CLI - Access chemical property data")
    
    parser.add_argument('--version', action='version', version=f'Chemical Data CLI v{__version__}')
    
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Search command
    search_parser = subparsers.add_parser('search', help='Search for chemicals by name')
    search_parser.add_argument('query', help='Chemical name to search for')
    search_parser.add_argument('--limit', type=int, default=10, help='Maximum number of results to show')
    
    # Lookup command
    lookup_parser = subparsers.add_parser('lookup', help='Look up properties for a specific chemical')
    lookup_parser.add_argument('identifier', help='Chemical name or CAS number')
    lookup_parser.add_argument('--all', action='store_true', help='Show all available properties')
    
    # Default to interactive mode if no args given
    if len(sys.argv) == 1:
        return None
        
    return parser.parse_args()

def main():
    """Main entry point for the application."""
    args = parse_args()
    
    if args is None:
        # Interactive mode
        print("Welcome to the Chemical Data CLI!")
        print("This tool provides access to chemical property data.")
        print("NOTE: For best results, use CAS numbers (e.g., 7732-18-5 for water)")
        main_menu()
    else:
        # Command line mode
        if args.command == 'search':
            search_chemicals(args.query, args.limit)
        elif args.command == 'lookup':
            if args.all:
                get_all_properties(args.identifier)
            else:
                lookup_chemical(args.identifier)

if __name__ == "__main__":
    main()