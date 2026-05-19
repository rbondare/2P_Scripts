import h5py
import numpy as np

def print_hdf5_item(name, obj):
    indent = name.count('/')
    padding = "  " * indent
    if isinstance(obj, h5py.Group):
        print(f"{padding}{name} (Group)")
    elif isinstance(obj, h5py.Dataset):
        print(f"{padding}{name} (Dataset): shape {obj.shape}, type {obj.dtype}")
        # If it's small, maybe show some info or if it's a list of refs
    
file_path = r'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat'

with h5py.File(file_path, 'r') as f:
    print("--- HDF5 Structure ---")
    f.visititems(print_hdf5_item)
    
    if 'CaData' in f:
        print("\n--- CaData Details ---")
        cadata = f['CaData']
        for key in cadata.keys():
            item = cadata[key]
            if isinstance(item, h5py.Dataset):
                print(f"{key}: {item.shape} {item.dtype}")
                if item.dtype == 'object': # Often references in v7.3
                    print(f"  (Likely references to other objects)")
            else:
                print(f"{key}: Group with {len(item.keys())} items")

