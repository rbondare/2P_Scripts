#!/usr/bin/env python3
"""
Interactive Axon Detection Parameter Tuner
Allows quick adjustment of detection parameters without editing code
"""

import numpy as np
import os
from pathlib import Path

def print_menu():
    print("\n" + "="*70)
    print("SUITE2P AXON DETECTION PARAMETER TUNER")
    print("="*70)
    print("\nAdjustment Presets:")
    print("  1. CONSERVATIVE  - Lower false positives, fewer detections")
    print("  2. BALANCED      - Medium sensitivity (default)")
    print("  3. AGGRESSIVE    - Higher sensitivity, more detections")
    print("  4. ULTRA-FINE    - For very faint, thin axons")
    print("  5. CUSTOM        - Enter your own parameters")
    print("  6. INFO          - Show current parameters")
    print("  0. EXIT          - Quit")
    print("\nChoice: ", end="")

def get_preset(choice):
    """Return parameter preset based on choice"""
    
    presets = {
        'CONSERVATIVE': {
            'cellprob_threshold': 0.2,
            'threshold_scaling': 0.6,
            'spatial_hp_detect': 12.0,
            'diameter': [2.5, 2.5],
            'block_size': [64.0, 64.0],
            'npix_norm_max': 80.0,
        },
        'BALANCED': {
            'cellprob_threshold': 0.1,
            'threshold_scaling': 0.5,
            'spatial_hp_detect': 10.0,
            'diameter': [3.0, 3.0],
            'block_size': [64.0, 64.0],
            'npix_norm_max': 100.0,
        },
        'AGGRESSIVE': {
            'cellprob_threshold': 0.05,
            'threshold_scaling': 0.35,
            'spatial_hp_detect': 8.0,
            'diameter': [2.5, 2.5],
            'block_size': [64.0, 64.0],
            'npix_norm_max': 120.0,
        },
        'ULTRA_FINE': {
            'cellprob_threshold': 0.02,
            'threshold_scaling': 0.25,
            'spatial_hp_detect': 5.0,
            'diameter': [2.0, 2.0],
            'block_size': [32.0, 32.0],
            'npix_norm_max': 60.0,
        },
    }
    
    return presets.get(choice, None)

def print_info(params):
    """Print current parameter settings"""
    print("\n" + "-"*70)
    print("CURRENT DETECTION PARAMETERS")
    print("-"*70)
    for key, val in sorted(params.items()):
        print(f"  {key:30s}: {val}")
    print("-"*70)

def custom_input():
    """Get custom parameters from user"""
    params = {}
    
    print("\nEnter custom parameters (or press Enter to skip):")
    print("-"*70)
    
    settings = [
        ('cellprob_threshold', 'Cell probability threshold (0.0-1.0, lower=more sensitive)', float),
        ('threshold_scaling', 'Threshold scaling (0.1-2.0, lower=more sensitive)', float),
        ('spatial_hp_detect', 'Spatial high-pass filter (5-25, lower=finer detail)', float),
        ('diameter_x', 'Diameter X (μm, typical: 2-5)', float),
        ('diameter_y', 'Diameter Y (μm, typical: 2-5)', float),
        ('block_size_x', 'Block size X (32-256, smaller=slower but finer)', float),
        ('block_size_y', 'Block size Y (32-256, smaller=slower but finer)', float),
        ('npix_norm_max', 'Max ROI pixels (30-150, smaller=thinner axons)', float),
    ]
    
    for key, desc, dtype in settings:
        user_input = input(f"\n  {desc}\n  [{key}]: ").strip()
        
        if user_input:
            try:
                val = dtype(user_input)
                if key == 'diameter_x':
                    params['diameter'] = [val, params.get('diameter', [3.0, 3.0])[1]]
                elif key == 'diameter_y':
                    params['diameter'] = [params.get('diameter', [3.0, 3.0])[0], val]
                elif key == 'block_size_x':
                    params['block_size'] = [val, params.get('block_size', [64.0, 64.0])[1]]
                elif key == 'block_size_y':
                    params['block_size'] = [params.get('block_size', [64.0, 64.0])[0], val]
                else:
                    params[key] = val
            except ValueError:
                print(f"  ✗ Invalid input, skipping {key}")
    
    return params

def apply_preset_to_file(preset_name, custom_params=None):
    """Apply preset to the axon detection script"""
    
    script_file = Path(r"C:\Users\rbondarenko\projects\2P_Scripts\suite2p_analysis\run_suite2p_RB_axons_V1.py")
    
    if not script_file.exists():
        print(f"✗ Script not found: {script_file}")
        return False
    
    # Get preset values
    if preset_name == 'CUSTOM' and custom_params:
        params = custom_params
    else:
        params = get_preset(preset_name)
    
    if not params:
        print("✗ Unknown preset")
        return False
    
    print("\n" + "-"*70)
    print("PARAMETERS TO APPLY:")
    print("-"*70)
    for key, val in sorted(params.items()):
        print(f"  {key:30s}: {val}")
    print("-"*70)
    
    confirm = input("\nApply these parameters? (y/n): ").strip().lower()
    if confirm != 'y':
        print("✗ Cancelled")
        return False
    
    # Read the script
    with open(script_file, 'r') as f:
        content = f.read()
    
    # Update parameters in the detection section
    updates = []
    
    if 'cellprob_threshold' in params:
        old = "'cellprob_threshold': 0.1,"
        new = f"'cellprob_threshold': {params['cellprob_threshold']},"
        content = content.replace(old, new)
        updates.append(f"cellprob_threshold: 0.1 → {params['cellprob_threshold']}")
    
    if 'threshold_scaling' in params:
        old = "'threshold_scaling': 0.5,"
        new = f"'threshold_scaling': {params['threshold_scaling']},"
        content = content.replace(old, new)
        updates.append(f"threshold_scaling: 0.5 → {params['threshold_scaling']}")
    
    if 'spatial_hp_detect' in params:
        old = "'spatial_hp_detect': 10.0,"
        new = f"'spatial_hp_detect': {params['spatial_hp_detect']},"
        content = content.replace(old, new)
        updates.append(f"spatial_hp_detect: 10.0 → {params['spatial_hp_detect']}")
    
    if 'diameter' in params:
        dia = params['diameter']
        old = "'diameter': [3.0, 3.0],"
        new = f"'diameter': [{dia[0]}, {dia[1]}],"
        content = content.replace(old, new)
        updates.append(f"diameter: [3.0, 3.0] → {dia}")
    
    if 'block_size' in params:
        block = params['block_size']
        old = "'block_size': [64.0, 64.0],"
        new = f"'block_size': [{block[0]}, {block[1]}],"
        content = content.replace(old, new)
        updates.append(f"block_size: [64.0, 64.0] → {block}")
    
    if 'npix_norm_max' in params:
        old = "'npix_norm_max': 100.0,"
        new = f"'npix_norm_max': {params['npix_norm_max']},"
        content = content.replace(old, new)
        updates.append(f"npix_norm_max: 100.0 → {params['npix_norm_max']}")
    
    # Write back the script
    with open(script_file, 'w') as f:
        f.write(content)
    
    print("\n✓ Script updated successfully!")
    print("\nChanges applied:")
    for update in updates:
        print(f"  - {update}")
    
    print(f"\nTo run with new parameters:")
    print(f"  python {script_file.name}")
    
    return True

if __name__ == "__main__":
    print("\n" + "="*70)
    print("SUITE2P AXON DETECTION PARAMETER TUNER")
    print("="*70)
    
    while True:
        print_menu()
        
        try:
            choice = input().strip()
            
            if choice == '0':
                print("\n✓ Goodbye!")
                break
            
            elif choice == '1':
                apply_preset_to_file('CONSERVATIVE')
            
            elif choice == '2':
                apply_preset_to_file('BALANCED')
            
            elif choice == '3':
                apply_preset_to_file('AGGRESSIVE')
            
            elif choice == '4':
                apply_preset_to_file('ULTRA_FINE')
            
            elif choice == '5':
                params = custom_input()
                apply_preset_to_file('CUSTOM', params)
            
            elif choice == '6':
                default_params = get_preset('BALANCED')
                print_info(default_params)
            
            else:
                print("✗ Invalid choice")
        
        except KeyboardInterrupt:
            print("\n\n✓ Interrupted")
            break
        
        except Exception as e:
            print(f"✗ Error: {e}")
