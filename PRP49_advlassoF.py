import os
import sys
import subprocess
import argparse
from collections import defaultdict
import random

def main():
    parser =argparse.ArgumentParser(description='为每个蛋白质随机选择random个lasso肽进行vina对接')
    parser.add_argument('--random', type=int, default=5, help='为每个蛋白进行随机抽样的个数')
    parser.add_argument('--wd', help='包含输入配体文件的目录')
    parser.add_argument('--rd', help='包含输入受体文件的目录')
    parser.add_argument('--script', default='PRP49_advlassoQ.py', help='对接流水线脚本路径')
    parser.add_argument('--ext', default='.pdb', help='输入文件的扩展名（默认为.pdb）')
    parser.add_argument('--od', help='输出目录（默认为输入目录）')
    parser.add_argument('--force', action='store_true', help='覆盖已存在的输出文件')

    args = parser.parse_args()

    # 检查输入目录
    wd = args.wd
    if not os.path.isdir(wd):
        print(f"错误：{wd} 不是有效目录")
        sys.exit(1)

    rd = args.rd
    if not os.path.isdir(rd):
        print(f"错误：{rd}不是有效目录")
        sys.exit(1)

    # 确定输出目录
    od = args.od if args.od else wd
    os.makedirs(od, exist_ok=True)

    # 获取所有指定扩展名的文件
    lps = [f for f in os.listdir(wd) if f.lower().endswith(args.ext.lower())]

    rps = defaultdict(list)

    rps_list = [f for f in os.listdir(rd) if f.lower().endswith(args.ext.lower())]

    sam_num = min(len(lps), args.random)
    for i in rps_list:
        if not sam_num:
            break
        selected_lps = random.sample(lps, sam_num)
        rps[i] = selected_lps

    for rp in rps:
        rppath = os.path.join(rd, rp)
        for lp in rps[rp]:
            base, _ = os.path.splitext(lp)
            rbase, _ = os.path.splitext(rp)
            oname = rbase + 'vs' + base +'.pdbqt'
            o = os.path.join(od, oname)

            # 检查输出文件
            if os.path.exists(o) and not args.force:
                print(f"跳过 {o}（输出文件已存在，使用 --force 覆盖）")
                continue

            print(f"处理{base}对接{rbase}...")
            cmd = [
                sys.executable,
                args.script,
                "--wd",wd,
                "--lp",lp,
                "--rp",rppath,
                "--o",o,
                "--compute"
                ]
            try:
                subprocess.run(cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"处理 {base}对接{rbase} 时出错：{e}")
                continue
    print("批量处理完成")

if __name__ == "__main__":
    main()

#python PRP49_advlassoF.py --wd C:\autodockvina\prp49\lassodata --rd C:\autodockvina\prp49\prodata --od C:\autodockvina\prp49\outdata