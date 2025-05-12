#!/usr/bin/env python3
"""
Validation utilities for chemical data processing.
"""
import re

def is_valid_cas(cas):
    """Check if a string appears to be a CAS number format"""
    pattern = r'^\d{1,7}-\d{2}-\d$'
    return bool(re.match(pattern, cas))