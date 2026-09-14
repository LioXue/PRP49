import os
import sys
import subprocess
import argparse
import re
import numpy as np

# 0.配置路径
parser = argparse.ArgumentParser(description='自动化分子对接流水线，必须输入--rp,--wd,--lp,--o')
# 路径参数（均提供默认值）
parser.add_argument('--mgltools-python', default=r"C:\Program Files (x86)\MGLTools-1.5.7\python.exe",
                    help='MGLTools Python 解释器路径')
parser.add_argument('--prepare-receptor', default=r"C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py",
                    help='prepare_receptor4.py 脚本路径')
parser.add_argument('--prepare-ligand', default=r"C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py",
                    help='prepare_ligand4.py 脚本路径')
parser.add_argument('--vina', default=r"C:\autodockvina\vina_1.2.7_win.exe",
                    help='AutoDock Vina 可执行文件路径')
parser.add_argument('--rp', required=True,
                    help='受体文件路径')
parser.add_argument('--rq',
                    help='受体MGLTool处理结果PDBQT文件路径')
parser.add_argument('--rforce', action='store_true',
                    help='覆盖受体已有PDBQT文件')
parser.add_argument('--wd', required=True,
                    help='工作目录（配体文件存放处）')
parser.add_argument('--lp', required=True,
                    help='配体PDB文件名（相对于工作目录）')
parser.add_argument('--ldq',default="default.pdbqt",
                    help='配体MGLTools处理结果PDBQT文件名（相对于工作目录）')
parser.add_argument('--notldforce',action='store_true',
                    help='不覆盖配体已有PDBQT文件')
parser.add_argument('--lq',
                    help='配体主链锁定PDBQT文件名（相对于工作目录），如无输入则使用lp的文件名命名，如启用notlock则无需输入')
parser.add_argument('--lforce',action='store_true',
                    help='覆盖配体已有PDBQT文件')
parser.add_argument('--c',default="config.txt",
                    help='Vina 配置文件文件名（相对于工作目录）')
parser.add_argument('--o', required=True,
                    help='对接结果PDBQT文件绝对路径')
parser.add_argument('--notlock', action='store_true',
                    help='停用 lasso 肽主链锁定')

# 对接盒子参数（可自定义）
parser.add_argument('--compute', action='store_true',
                    help='启用计算')
parser.add_argument('--x', type=float, default=0,
                    help='盒子中心 X 坐标')
parser.add_argument('--y', type=float, default=0,
                    help='盒子中心 Y 坐标')
parser.add_argument('--z', type=float, default=-0,
                    help='盒子中心 Z 坐标')
parser.add_argument('--size-x', type=float, default=30.0,
                    help='盒子 X 方向尺寸')
parser.add_argument('--size-y', type=float, default=30.0,
                    help='盒子 Y 方向尺寸')
parser.add_argument('--size-z', type=float, default=30.0,
                    help='盒子 Z 方向尺寸')
parser.add_argument('--exhaustiveness', type=int, default=3,
                    help='计算精度（默认3）')
parser.add_argument('--num-modes', type=int, default=1,
                    help='输出构象数（默认1）')
parser.add_argument('--energy-range', type=float, default=1,
                    help='能量范围（默认1）')

args = parser.parse_args()

WORK_DIR = args.wd

LIGAND_PDB = args.lp   # 相对路径
LIGAND_DEFAULT_PDBQT = args.ldq   # 相对路径
if not args.lq:
    lq = os.path.splitext(args.lp)[0] + ".pdbqt"
    LIGAND_PDBQT = os.path.join(WORK_DIR, lq)
else:
    LIGAND_PDBQT = os.path.join(WORK_DIR, args.lq)
CONFIG_FILE = os.path.join(WORK_DIR, args.c)

if not args.rq:
    rq = os.path.splitext(args.rp)[0] + ".pdbqt"
else:
    rq = args.rq

# 1.受体预处理
print("MGLTools处理受体...")
if os.path.exists(rq) and not args.rforce:
    print(f"受体文件已存在，跳过受体处理")
else:
    print("生成受体文件...")
    cmd = [
        args.mgltools_python,
        args.prepare_receptor,
        "-r", args.rp,
        "-o", rq,
        "-A", "hydrogens",
        "-U", "waters"
    ]
    subprocess.run(cmd, check=True)
    print("受体文件生成完成")

# 2.配体预处理
print("MGLTools处理配体...")
if os.path.exists(os.path.join(WORK_DIR,LIGAND_DEFAULT_PDBQT)) and args.notldforce:
    print(f"配体文件已存在: {LIGAND_DEFAULT_PDBQT}，跳过配体处理。")
else:
    cmd = [
        args.mgltools_python,
        args.prepare_ligand,
        "-l", LIGAND_PDB,
        "-o", LIGAND_DEFAULT_PDBQT,
        "-A", "bonds_hydrogens",
        "-R", "1"
        ]
    subprocess.run(cmd, cwd=WORK_DIR, check=True)
    print("配体文件生成完成。")

LIGAND_DEFAULT_PDBQT_full = os.path.join(WORK_DIR, LIGAND_DEFAULT_PDBQT)

# 3.配体主链锁定
if not args.notlock:
    if os.path.exists(LIGAND_PDBQT) and not args.lforce:
        print(f"配体文件已存在: {LIGAND_PDBQT}，跳过配体锁定。")
    else:
        print("锁定主链扭转键...")
        trans_script_dir = os.path.dirname(__file__)
        trans_script = os.path.join(trans_script_dir, "PRP49_translasso.py")
        cmd = [
            sys.executable,
            trans_script,
            LIGAND_DEFAULT_PDBQT_full,
            LIGAND_PDBQT
            ]
        subprocess.run(cmd, check=True)
        print("配体锁定完成。")
else:
    LIGAND_PDBQT = LIGAND_DEFAULT_PDBQT_full


# 4.配置vina config
print("创建 Vina 配置文件...")

def read_position(info):
    """读取原子位置"""
    pattern = r'^ATOM\s+.*?\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+).*\s([A-Za-z]+)\s*$'
    match = re.match(pattern, info)
    if match:
        x, y, z, atom = match.groups()
        return float(x), float(y), float(z), atom
    return None

if args.compute:
    v_list=[]
    with open(rq, "r") as f:
        for line in f:
            if read_position(line):
                x, y, z, atom = read_position(line)
                x_max, y_max, z_max = x, y, z
                x_min, y_min, z_min = x, y, z
                break
            else:
                continue
        for line in f:
            if not read_position(line):
                continue
            xs, ys, zs, atom = read_position(line)
            x, y, z = xs, ys, zs
            if x>x_max:
                x_max = x
            if y>y_max:
                y_max = y
            if z>z_max:
                z_max = z
            if x<x_min:
                x_min = x
            if y<y_min:
                y_min = y
            if z<z_min:
                z_min = z
            if atom == "C":
                v = [x,y,z]
                v_list.append(v)
    center = np.array(v_list).mean(axis=0) #重心机制，center位于蛋白质的重心
    xc = center[0]
    yc = center[1]
    zc = center[2]
    delta_x = abs(max(xc-x_min, xc-x_max, key=lambda t: abs(t))) #距离重心更远的一边会获得更大的空间
    delta_y = abs(max(yc-y_min, yc-y_max, key=lambda t: abs(t)))
    delta_z = abs(max(zc-z_min, zc-z_max, key=lambda t: abs(t)))
    xsize = 2*delta_x
    ysize = 2*delta_y
    zsize = 2*delta_z
else:
    xc = args.x
    yc = args.y
    zc = args.z
    xsize = args.size_x
    ysize = args.size_y
    zsize = args.size_z
if args.compute and xsize*ysize*zsize > 125000:
    exh = max(5,args.exhaustiveness)
else:
    exh = args.exhaustiveness
config_content = f"""receptor = {rq}
ligand = {LIGAND_PDBQT}
out = {args.o}
center_x = {xc}
center_y = {yc}
center_z = {zc}
size_x = {xsize}
size_y = {ysize}
size_z = {zsize}
exhaustiveness = {exh}
num_modes = {args.num_modes}
energy_range = {args.energy_range}
"""
with open(CONFIG_FILE, 'w') as f:
    f.write(config_content)
print("配置文件已创建。")

# 5.运行vina
print("运行vina...")
cmd = [args.vina, "--config", CONFIG_FILE]
subprocess.run(cmd, check=True)
print("对接已完成")