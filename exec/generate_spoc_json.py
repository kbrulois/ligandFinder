import os
import glob
import pickle
import argparse
import json
import numpy as np

def convert_pkl_to_json(pkl_file: str, mode: str) -> None:
    """Load an AlphaFold2 result_model_*.pkl file and save pae/ptm/iptm to a .json."""
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)
    
    required_keys = ["predicted_aligned_error", "ptm", "iptm"]
    if not all(k in data for k in required_keys):
        print(f"Skipping {pkl_file} - missing required keys.")
        return

    # Convert NumPy arrays/floats into vanilla Python structures
    json_dict = {
        "pae": data["predicted_aligned_error"].tolist(),
        "ptm": float(data["ptm"]),
        "iptm": float(data["iptm"])
    }
    
    if mode == "og":
      json_file = os.path.splitext(pkl_file)[0].replace("result_", "paeSPOC_") + ".json"
    if mode == "rn":
      json_file = os.path.splitext(pkl_file)[0].replace("_res_", "_spc_") + ".json"
    with open(json_file, "w") as f:
        json.dump(json_dict, f)
    print(f"Created {json_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .pkl files to .json with pae/ptm/iptm")
    parser.add_argument("input", type=str, help="Path to run directory containing AF output directories")
    parser.add_argument("mode", type=str, help="original or renamed")
    args = parser.parse_args()

    pattern = os.path.join(args.input, "*/*.pkl")
    pkl_files = glob.glob(pattern)

    for pkl_file in pkl_files:
        convert_pkl_to_json(pkl_file=pkl_file, mode=args.mode)
