import os
import re
import matplotlib.pyplot as plt
import numpy as np

base_dir = "."
colors = ['#FFFFFF', '#D3D3D3', '#808080', '#000000', '#A9A9A9', '#696969']

def parse_csv_for_error(filepath):
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('4-cycles,'):
                    parts = line.strip().split(',')
                    return float(parts[1])
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
    return None

def main():
    for dataset_name in os.listdir(base_dir):
        dataset_path = os.path.join(base_dir, dataset_name)
        if not os.path.isdir(dataset_path) or dataset_name.startswith('.'):
            continue

        data = {}
        epsilons_set = set()
        

        for filename in os.listdir(dataset_path):
            if not filename.endswith('.csv'):
                continue

            match = re.search(r'alg(\d+)_eps([\d\.]+)', filename)
            if match:
                alg = f"alg{match.group(1)}"
                eps = float(match.group(2))
                
                filepath = os.path.join(dataset_path, filename)
                error_val = parse_csv_for_error(filepath)
                
                if error_val is not None:
                    if alg not in data:
                        data[alg] = {}
                    data[alg][eps] = error_val
                    epsilons_set.add(eps)

        if not data:
            continue
            
        print(f"Plotting: {dataset_name} ...")

        sorted_eps = sorted(list(epsilons_set))
        sorted_algs = sorted(list(data.keys()))
        
        x = np.arange(len(sorted_eps))
        num_algs = len(sorted_algs)
        total_width = 0.8
        width = total_width / num_algs
        fig, ax = plt.subplots(figsize=(10, 6))

        for i, alg in enumerate(sorted_algs):
            y_values = [data[alg].get(eps, np.nan) for eps in sorted_eps]
            
            offset = x + (i * width) - (total_width / 2) + (width / 2)
            
            ax.bar(offset, y_values, width, label=alg, 
                   color=colors[i % len(colors)], edgecolor='black', zorder=3)

        ax.set_yscale('log')
        ax.set_xlabel('Epsilon (\u03B5)', fontsize=12)
        ax.set_ylabel('Mean Relative Error', fontsize=12)
        ax.set_title(f'Performance on {dataset_name} Dataset', fontsize=14)
        
        ax.set_xticks(x)
        ax.set_xticklabels([str(eps) for eps in sorted_eps])
        
        ax.legend(loc='upper left', frameon=False, ncol=num_algs)
        
        plt.tight_layout()
        output_filename = f"plot_{dataset_name}.png"
        plt.savefig(output_filename, dpi=300)
        plt.close()
        print(f"Generated: {output_filename}")

if __name__ == "__main__":
    main()