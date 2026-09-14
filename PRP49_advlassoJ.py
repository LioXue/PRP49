import os
import sys
import subprocess
import argparse
import random

def main():
    parser =argparse.ArgumentParser(description='为每个蛋白质随机选择lasso肽vina成功对接random次')
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

    ext = args.ext if args.ext.startswith('.') else '.' + args.ext

    lps = [f
           for f in os.listdir(wd)
           if os.path.isfile(os.path.join(wd, f)) and f.lower().endswith(ext.lower())]

    rps_list = [f for f in os.listdir(rd) if f.lower().endswith(ext.lower())]

    sam_num = max(0,min(len(lps), args.random))

    for rp in rps_list:
        rppath = os.path.join(rd, rp)
        lps_list = lps[:]  # 每个 rp 都重新洗牌，独立抽样
        random.shuffle(lps_list)

        idx = 0  # 已经获得的结果数
        while idx < sam_num and lps_list:
            lp = lps_list.pop(0)
            base, _ = os.path.splitext(lp)
            rbase, _ = os.path.splitext(rp)
            oname = rbase + 'vs' + base + '.pdbqt'
            o = os.path.join(od, oname)

            # 情况1：输出文件已存在且不覆盖 —— 它本身就是一个结果
            if os.path.exists(o) and not args.force:
                print(f"跳过 {o}（输出文件已存在，使用 --force 覆盖）")
                idx += 1
                continue

            # 情况2：正常对接
            print(f"处理 {base} 对接 {rbase} ...")
            cmd = [
                sys.executable,
                args.script,
                "--wd", wd,
                "--lp", lp,
                "--rp", rppath,
                "--o", o,
                "--compute"
            ]
            try:
                subprocess.run(cmd, check=True)
                idx += 1  # 成功，结果数 +1
            except subprocess.CalledProcessError as e:
                print(f"处理 {base} 对接 {rbase} 时出错：{e}")
                # 失败不增加 idx，继续尝试下一个配体
                continue

        # 配体用尽但结果数不够，给出警告
        if idx < sam_num:
            print(f"警告：{rp} 只获得 {idx}/{sam_num} 个结果（配体已用尽或全部失败）")
    print("批量处理完成")

if __name__ == "__main__":
    main()

