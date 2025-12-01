from io import TextIOWrapper
from graph_generator import *
from typing import Type
import os
from datetime import datetime
import shutil
import tqdm

def run(generator: Type[GraphGenerator]):
    print(f"Running experiment for generator: {generator.name}")
    t_saved = [1, 10, 100, 1000]
    dir_path = f"results/{generator.dir_name}"
    file_template = dir_path + "/" + "timeseries_t{}.txt"
    os.makedirs(dir_path, exist_ok=True)

    for t in t_saved:
        file_path = file_template.format(t)
        if os.path.exists(file_path):
            print(f"WARNING: Ja existeix el fitxer {file_path}.")
            old_dir = os.path.join(dir_path, "old")
            os.makedirs(old_dir, exist_ok=True)
            mtime = os.path.getmtime(file_path)
            timestamp = datetime.fromtimestamp(mtime).strftime("%Y%m%d-%H%M%S")
            base = os.path.splitext(os.path.basename(file_path))[0]
            new_path = os.path.join(old_dir, f"{base}-{timestamp}.txt")
            shutil.move(file_path, new_path)
            print(f"Moved existing file to {new_path}.")
    files: dict[int, TextIOWrapper] = {}
    try:
        for t in t_saved:
            file_path = file_template.format(t)
            files[t] = open(file_path, "w")

        for T, G in tqdm.tqdm(generator.generate(), total=Constants.T_MAX):
            for t, f in files.items():
                if G.has_node(t):
                    f.write(f"{T} {G.degree[t]}\n")
        last_G = G
    finally:
        for f in files.values():
            f.close()
            print(f"Closed file {f.name}.")

    G = last_G
    degrees = sorted([d for n, d in G.degree()], reverse=True)
    with open(os.path.join(dir_path, "degrees_tmax.txt"), "w") as f:
        for d in degrees:
            f.write(f"{d}\n")
        print(f"Saved degree distribution to {f.name}.")

if __name__ == "__main__":
    classes = [BarabasiAlbert, BA_RandomAttachment, BA_StaticNodes]
    for cls in classes:
        run(cls)