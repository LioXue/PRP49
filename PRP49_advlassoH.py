import os
import sys
import subprocess
import argparse
import random

def main():
    parser =argparse.ArgumentParser(description='为目标蛋白质随机选择random个lasso肽进行vina对接')
    parser.add_argument('--random', type=int, default=5, help='为每个蛋白进行随机抽样的个数')
    parser.add_argument('--wd', help='包含输入配体文件的目录')
    parser.add_argument('--rp', help='受体文件路径')
    parser.add_argument('--rq', help='受体PDBQT文件路径')
    parser.add_argument('--ext', default='.pdb', help='输入文件扩展名（默认 .pdb）')
    parser.add_argument('--script', default='PRP49_advlassoQ.py', help='对接流水线脚本路径')
    parser.add_argument('--od', help='输出目录（默认为输入目录）')
    parser.add_argument('--force', action='store_true', help='覆盖已存在的输出文件')

    args = parser.parse_args()

    # 检查输入目录
    wd = args.wd
    if not os.path.isdir(wd):
        print(f"错误：{wd} 不是有效目录")
        sys.exit(1)

    # 确定受体
    rp = args.rp
    rq = args.rq
    rname = os.path.splitext(os.path.basename(rp))[0]

    # 确定输出目录
    od = args.od if args.od else wd
    os.makedirs(od, exist_ok=True)

    # 获取所有指定扩展名的文件
    ext = args.ext if args.ext.startswith('.') else '.' + args.ext
    ext = ext.lower()

    lps = [
        f for f in os.listdir(wd)
        if os.path.isfile(os.path.join(wd, f)) and f.lower().endswith(ext.lower())
    ]

    sam_num = max(0,min(len(lps),args.random))

    lps = random.sample(lps,sam_num)

    for lp in lps:
        base, _ = os.path.splitext(lp)
        lq = base + '.pdbqt'
        dq = f"{rname}vs{base}" + '.pdbqt'
        o = os.path.join(od, dq)

        # 跳过已存在且未强制覆盖的文件
        if os.path.exists(o) and not args.force:
            print(f"跳过 {o}（输出文件已存在，使用 --force 覆盖）")
            continue

        print(f"处理 {lp} -> {o}...")
        # 调用对接流水线脚本
        cmd = [
            sys.executable,
            args.script,
            "--wd",wd,
            "--lp",lp,
            "--lq",lq,
            "--o",o,
            "--rp",rp,
            "--rq",rq,
            "--compute"
        ]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"处理 {lp} 时出错：{e}")
            continue
    print("批量处理完成")

if __name__ == "__main__":
    main()