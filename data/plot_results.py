import os
import re
import matplotlib.pyplot as plt
import numpy as np

# 设置工作目录为脚本所在目录 (即 ~/Triangle4CycleShuffle/data)
base_dir = "."

# 图三风格的灰度配色方案
colors = ['#FFFFFF', '#D3D3D3', '#808080', '#000000', '#A9A9A9', '#696969']

def parse_csv_for_error(filepath):
    """读取 CSV 文件并提取 4-cycles 的 AVG(rel-err)"""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # 寻找目标行
                if line.startswith('4-cycles,'):
                    parts = line.strip().split(',')
                    # 返回第二个元素，即 AVG(rel-err)
                    return float(parts[1])
    except Exception as e:
        print(f"读取文件出错 {filepath}: {e}")
    return None

def main():
    # 遍历当前目录下的所有项
    for dataset_name in os.listdir(base_dir):
        dataset_path = os.path.join(base_dir, dataset_name)
        
        # 只处理文件夹（即各个数据集，如 musae, Gplus 等）
        if not os.path.isdir(dataset_path) or dataset_name.startswith('.'):
            continue

        # 数据结构：data[alg_name][epsilon_value] = error_value
        data = {}
        epsilons_set = set()
        
        # 遍历数据集文件夹内的所有 CSV 文件
        for filename in os.listdir(dataset_path):
            if not filename.endswith('.csv'):
                continue
                
            # 使用正则表达式提取 alg 和 epsilon
            # 文件名示例: res_n22470_alg1_eps0.5-8_pair-1_itr20-1.csv
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

        # 如果这个文件夹里没有提取到有效数据，跳过
        if not data:
            continue
            
        print(f"正在绘制数据集: {dataset_name} ...")

        # 准备绘图数据，保证 x 轴 (epsilon) 和 图例 (算法) 是有序的
        sorted_eps = sorted(list(epsilons_set))
        sorted_algs = sorted(list(data.keys()))
        
        x = np.arange(len(sorted_eps))  # x 轴基准位置
        num_algs = len(sorted_algs)
        total_width = 0.8  # 柱状图总宽度
        width = total_width / num_algs  # 每个算法柱子的宽度
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 遍历算法画柱状图
        for i, alg in enumerate(sorted_algs):
            # 获取该算法在每个 epsilon 下的 error，如果没有数据则用 NaN 填充占位
            y_values = [data[alg].get(eps, np.nan) for eps in sorted_eps]
            
            # 计算每个柱子的偏移量，使其居中对称排列
            offset = x + (i * width) - (total_width / 2) + (width / 2)
            
            # 绘制柱状图，使用灰度配色，加上黑色边框以模仿图三风格
            ax.bar(offset, y_values, width, label=alg, 
                   color=colors[i % len(colors)], edgecolor='black', zorder=3)

        # 设置图表格式 (对数坐标轴，模仿图三)
        ax.set_yscale('log')
        ax.set_xlabel('Epsilon (\u03B5)', fontsize=12)
        ax.set_ylabel('Mean Relative Error', fontsize=12)
        ax.set_title(f'Performance on {dataset_name} Dataset', fontsize=14)
        
        # 设置 x 轴刻度标签
        ax.set_xticks(x)
        ax.set_xticklabels([str(eps) for eps in sorted_eps])
        
        # 添加图例和网格线
        ax.legend(loc='upper left', frameon=False, ncol=num_algs) # 图例在一排显示
        
        # 调整布局并保存图片
        plt.tight_layout()
        output_filename = f"plot_{dataset_name}.png"
        plt.savefig(output_filename, dpi=300)
        plt.close()
        print(f"已生成图像: {output_filename}")

if __name__ == "__main__":
    main()