import scipy.io
import numpy as np

def print_struct(obj, name="", indent=0):
    tab = "  " * indent
    if isinstance(obj, np.ndarray):
        if obj.dtype.names is not None:
            print(f"{tab}{name} (struct array): shape {obj.shape}")
            for n in obj.dtype.names:
                # Inspecting first element for variety if it's an array of structs
                if obj.size > 0:
                    print_struct(obj[0][n], n, indent + 1)
                else:
                    print(f"{tab}  {n} (empty)")
        else:
            print(f"{tab}{name} (ndarray): shape {obj.shape}, dtype {obj.dtype}")
    elif isinstance(obj, dict):
        print(f"{tab}{name} (dict):")
        for k, v in obj.items():
            if not k.startswith('__'):
                print_struct(v, k, indent + 1)
    else:
        print(f"{tab}{name}: {type(obj)}")

file_path = r'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat'
data = scipy.io.loadmat(file_path, struct_as_record=False, squeeze_me=True)

print("--- File Contents ---")
for key in data:
    if not key.startswith('__'):
        val = data[key]
        if hasattr(val, '_fieldnames'):
            print(f"{key} (struct): fields {val._fieldnames}")
            for field in val._fieldnames:
                field_val = getattr(val, field)
                if isinstance(field_val, np.ndarray):
                     print(f"  {field}: shape {field_val.shape}")
                else:
                     print(f"  {field}: {type(field_val)}")
        elif isinstance(val, np.ndarray):
            print(f"{key}: array shape {val.shape}")
        else:
            print(f"{key}: {type(val)}")

if 'CaData' in data:
    print("\n--- CaData Structure ---")
    cadata = data['CaData']
    if hasattr(cadata, '_fieldnames'):
        for field in cadata._fieldnames:
             field_val = getattr(cadata, field)
             if isinstance(field_val, np.ndarray):
                  print(f"{field}: shape {field_val.shape}")
             else:
                  print(f"{field}: type {type(field_val)}")
    else:
        print(f"CaData is type {type(cadata)}")

