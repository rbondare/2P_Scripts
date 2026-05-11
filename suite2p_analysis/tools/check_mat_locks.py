import os, glob
plane_dir = r"Z:\group\joeschgrp\Group Members\Rima\DATA_2P\AnimalRB19_260320_1444\suite2p"
for f in glob.glob(os.path.join(plane_dir, "*.npy")):
    try:
        os.rename(f, f + ".testlock")
        os.rename(f + ".testlock", f)
        print("writable:", f)
    except PermissionError:
        print("locked (PermissionError):", f)
    except Exception as e:
        print("error:", f, type(e).__name__, e)