import re
import math
import PRP49_amino_acid_sidechain as aas
import PRP49_TailBackbone as tb
import argparse
import os

par = argparse.ArgumentParser(description="Translate a default PDBQT file to a lasso-friendly PDBQT")
par.add_argument("input_pdbqt", help="default PDBQT file to translate")
par.add_argument("output_pdbqt", help="translated PDBQT file")
args = par.parse_args()

input_pdbqt = args.input_pdbqt
output_pdbqt = args.output_pdbqt

#测试代码
#WORK_DIR = r"C:\autodockvina\prp49"
#input_pdbqt = os.path.join(WORK_DIR,"2LTI_test7.pdbqt") # 配体中间pdbqt(只要相对路径)
#output_pdbqt = os.path.join(WORK_DIR, "2LTI_test11.pdbqt") # 配体最终pdbqt

class atom_info:
    """储存原子信息"""
    def __init__(self, atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                 n=False, ca = False, c = False, o = False, h = False, oxt = False,
                 root = False, tail = False, idx = 0):
        self.atomnum1 = atomnum1
        self.atom = atom
        self.atomtip = atomtip
        self.atomnum2 = atomnum2
        self.residue = residue
        self.chain = chain
        self.num = num
        self.info = info
        self.gap = gap
        self.n = n
        self.ca = ca
        self.c = c
        self.o = o
        self.h = h
        self.oxt = oxt
        self.root = root
        self.tail = tail
        self.idx = idx
        self.sidechain = []

    def make_position(self):
        return position(self)

class position:
    """为原子信息节点建立位置"""
    def __init__(self, atom_info):
        self.atom_info = atom_info


def parse(line):
    pattern = r'^ATOM\s+\d+(\s+)(\d*)([A-Za-z])([A-Za-z]*)(\d*)(\s+)(\w+)(\s+)(\w*)(\s+)(\d+)(\s+)(.*)$'
    #识别ATOM     76 1HE2 GLN A  22       0.375   2.036   8.050  1.00  0.00     0.159 HD
    match = re.match(pattern, line)
    if match:
        gap1, atomnum1, atom, atomtip, atomnum2, gap2, residue, gap3, chain, gap4, num, gap5, info = match.groups()
        return atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, [gap1, gap2, gap3, gap4, gap5]
    return None


atom_info_list = [] #储存原子节点

with open(input_pdbqt,"r") as f:
    for line in f:
        parsed = parse(line)
        if parsed:
            atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap = parsed
            if atom == "N" and atomtip in ("", None):
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                                                n=True))
            elif atom == "C" and atomtip == "A":
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                                                ca=True))
            elif atom == "C" and atomtip in ("", None):
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                                                c=True))
            elif atom == "O" and atomtip in ("", None):
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                                                o=True))
            elif atom == "H" and atomtip in ("","N",None) and atomnum2 in ("","1",None):
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                                                h=True))
            elif atom == "O" and atomtip == "XT":
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap,
                                                oxt=True))
            else :
                atom_info_list.append(atom_info(atomnum1, atom, atomtip, atomnum2, residue, chain, num, info, gap))

residue_info_list = [] #list储存dict

for atom_info in atom_info_list:
    if atom_info.ca:
        residue_info_dict = {"ca":atom_info.make_position(),
                             "h":None, "n":None, "o":None, "c":None, "cb":None, "oxt":None,
                             "cg":[], "cd":[], "ce":[], "cz":[],"ox":[],"nx":[],"hx":[],"s":[]}
        # dict储存{原子：原子节点的位置，原子：位置}
        for i in atom_info_list:
            if i.num == atom_info.num and i.residue == atom_info.residue:
                if i.h:
                    residue_info_dict["h"] = i.make_position()
                elif i.n:
                    residue_info_dict["n"] = i.make_position()
                elif i.o:
                    residue_info_dict["o"] = i.make_position()
                elif i.c:
                    residue_info_dict["c"] = i.make_position()
                elif i.oxt:
                    residue_info_dict["oxt"] = i.make_position()
                elif i.atom == "C" and i.atomtip == "B":
                    residue_info_dict["cb"] = i.make_position()
                elif i.atom == "C" and i.atomtip == "G":
                    residue_info_dict["cg"].append(i.make_position())
                elif i.atom == "C" and i.atomtip == "D":
                    residue_info_dict["cd"].append(i.make_position())
                elif i.atom == "C" and i.atomtip == "E":
                    residue_info_dict["ce"].append(i.make_position())
                elif i.atom == "C" and i.atomtip == "Z":
                    residue_info_dict["cz"].append(i.make_position())
                elif i.atom == "O" and not i.o and not i.oxt:
                    residue_info_dict["ox"].append(i.make_position())
                elif i.atom == "N" and not i.n:
                    residue_info_dict["nx"].append(i.make_position())
                elif i.atom == "H" and not i.h:
                    residue_info_dict["hx"].append(i.make_position())
                elif i.atom == "S":
                    residue_info_dict["s"].append(i.make_position())
        aas.sidechain_find(residue_info_dict)
        residue_info_list.append(residue_info_dict) #将dict加进list

residue_info_list.sort(key=lambda item: int(item["ca"].atom_info.num)) #根据残基位数排序

def read_position(info):
    """读取原子位置"""
    pattern = r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+.*$"
    #识别0.375   2.036   8.050  1.00  0.00     0.159 HD
    match = re.match(pattern, info)
    if match:
        x, y, z = match.groups()
        return x, y, z
    return None

#找到环结束的位置
first_position = read_position(residue_info_list[0]["n"].atom_info.info)
v1 = [float(i) for i in first_position]
min_distance = 1000
ring_end_index = 0
if first_position:
    for j,i in enumerate(residue_info_list):
        if i["ca"].atom_info.residue in ("ASP","ASN","ASX","GLU","GLN","GLX") and j <= 0.6*len(residue_info_list):
            second_position = read_position(i["cd"][0].atom_info.info) if i["cd"] else read_position(i["cg"][0].atom_info.info)
            v2 = [float(i) for i in second_position]
            distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v2)))
            if distance < min_distance:
                min_distance = distance
                ring_end_index = j
            else:
                continue
        else:
            continue


#ring_end_index = 8

#找到尾巴开始的地方
tail_start_index = ring_end_index + 1
min_distance = 100000
ring_atom_position = [read_position(residue_info_list[i]["ca"].atom_info.info) for i in range(ring_end_index+1)]
other_atom_position = [read_position(residue_info_list[i]["ca"].atom_info.info) for i in range(ring_end_index+1, len(residue_info_list))]
for l, i in enumerate(other_atom_position):
    v = [float(j) for j in i]
    distance = 0
    for k in ring_atom_position:
        v3 = [float(j) for j in k]
        distance += math.sqrt(sum((a - b) ** 2 for a, b in zip(v, v3)))
    if distance < min_distance:
        min_distance = distance
        tail_start_index = ring_end_index + l + 1
    else:
        continue

print("calculation done, the iso-peptide position is:", ring_end_index+1, "the tail start position is:", tail_start_index+1)

#为tail_start索引前的所有主链原子打上root，为tail_start索引的N,H原子打上root
for i in range(tail_start_index + 1):
    for j in residue_info_list[i]:
        if residue_info_list[i][j] and j in ("n","h","ca","c","o") and i < tail_start_index:
            residue_info_list[i][j].atom_info.root = True
        elif residue_info_list[i][j] and j in ("n","h") and i == tail_start_index:
            residue_info_list[i][j].atom_info.root = True
        else:
            continue

for j in residue_info_list[ring_end_index]:
    if residue_info_list[ring_end_index][j] and j == "cb":
        residue_info_list[ring_end_index][j].atom_info.root = True
    elif residue_info_list[ring_end_index][j] and j in ("cg", "cd", "nx", "ox", "hx"):
        for x in residue_info_list[ring_end_index][j]:
            x.atom_info.root = True

#为tail_starts索引的CA,C,O及所有后续主链原子打上tail
for i in range(tail_start_index, len(residue_info_list)):
    for j in residue_info_list[i]:
        if residue_info_list[i][j] and j in ("ca","c","o") and i == tail_start_index:
            residue_info_list[i][j].atom_info.tail = True
        elif residue_info_list[i][j] and j in ("n","h","ca","c","o") and i > tail_start_index:
            residue_info_list[i][j].atom_info.tail = True
        else:
            continue

root = []

#将所有写有root的原子加入root_list
index = 1
branch_amino_acid = 0
branch_index = 0 #为tail起始提供索引
for i in atom_info_list:
    if i.root:
        root.append("ATOM"+" "*max(0,(7-len(str(index))))+f"{index}{i.gap[0]}{i.atomnum1}{i.atom}{i.atomtip}{i.atomnum2}{i.gap[1]}{i.residue}{i.gap[2]}{i.chain}{i.gap[3]}{i.num}{i.gap[4]}{i.info}\n")
        i.idx = index
        if i.atom == "N" and int(i.num) > branch_amino_acid:
            branch_amino_acid = int(i.num)
            branch_index = i.idx
        index +=1
    else:
        continue

tail = tb.TailBackbone()#建立TB数据结构

for j in residue_info_list:
    if j["ca"].atom_info.tail:
        tail.add(j)
    else:
        continue

def write_sidechain(base_index, position, file, index):
    new_index = index
    if position.atom_info.sidechain:
        if not isinstance(position.atom_info.sidechain, list) and not position.atom_info.sidechain.atom_info.root:
            file.write("BRANCH" + " " * max(0, (4 - len(str(base_index)))) + f"{base_index}" + " " * max(0, (4 - len(str(base_index)))) + f"{new_index}\n")
            record_tuple = base_index, new_index
            i = position.atom_info.sidechain.atom_info
            file.write("ATOM" + " " *max(0,(7-len(str(new_index))))+f"{new_index}{i.gap[0]}{i.atomnum1}{i.atom}{i.atomtip}{i.atomnum2}{i.gap[1]}{i.residue}{i.gap[2]}{i.chain}{i.gap[3]}{i.num}{i.gap[4]}{i.info}\n")
            i.idx = new_index
            new_index += 1
            if i.sidechain:
                new_index = write_sidechain(i.idx, i, file, new_index)
            else:
                pass
            file.write("ENDBRANCH"+" "*max(0,(4-len(str(record_tuple[0]))))+f"{record_tuple[0]}"+" "*max(0, (4 - len(str(record_tuple[1]))))+f"{record_tuple[1]}\n")
        elif isinstance(position.atom_info.sidechain, list) and not position.atom_info.sidechain[0].atom_info.root:
            file.write("BRANCH" + " " * max(0, (4 - len(str(base_index)))) + f"{base_index}" + " " * max(0, (4 - len(str(base_index)))) + f"{new_index}\n")
            record_tuple = base_index, new_index
            for j in position.atom_info.sidechain:
                i = j.atom_info
                file.write("ATOM" + " " *max(0,(7-len(str(new_index))))+f"{new_index}{i.gap[0]}{i.atomnum1}{i.atom}{i.atomtip}{i.atomnum2}{i.gap[1]}{i.residue}{i.gap[2]}{i.chain}{i.gap[3]}{i.num}{i.gap[4]}{i.info}\n")
                i.idx = new_index
                new_index += 1
            for j in position.atom_info.sidechain:
                if j.atom_info.sidechain:
                    new_index = write_sidechain(j.atom_info.idx, j, file, new_index)
                else:
                    pass
            file.write("ENDBRANCH"+" "*max(0,(4-len(str(record_tuple[0]))))+f"{record_tuple[0]}"+" "*max(0, (4 - len(str(record_tuple[1]))))+f"{record_tuple[1]}\n")
        else:
            pass
        return new_index
    else:
        return new_index


def write_tail(branch_index, tail, file, index):
    base_index = branch_index
    new_index = index
    endbranch_list = []
    for j in tail:
        if j:
            if not isinstance(j, tuple):
                i = j.atom_info
                file.write("BRANCH" + " " * max(0, (4 - len(str(base_index)))) + f"{base_index}"+" "*max(0, (4 - len(str(base_index))))+f"{new_index}\n")
                endbranch_list.append((base_index, new_index))
                file.write("ATOM" + " " *max(0,(7-len(str(new_index))))+f"{new_index}{i.gap[0]}{i.atomnum1}{i.atom}{i.atomtip}{i.atomnum2}{i.gap[1]}{i.residue}{i.gap[2]}{i.chain}{i.gap[3]}{i.num}{i.gap[4]}{i.info}\n")
                i.idx = new_index
                base_index = new_index
                new_index += 1
                new_index = write_sidechain(base_index, j, file, new_index)
            elif len(j) == 4:
                c, o, n, h = j
                conh = [c.atom_info, o.atom_info, n.atom_info, h.atom_info]
                file.write("BRANCH" + " " * max(0,(4 - len(str(base_index)))) + f"{base_index}"+" "*max(0, (4 - len(str(base_index))))+f"{new_index}\n")
                endbranch_list.append((base_index, new_index))
                for x in conh:
                    file.write("ATOM"+ " "*max(0,(7-len(str(new_index))))+f"{new_index}{x.gap[0]}{x.atomnum1}{x.atom}{x.atomtip}{x.atomnum2}{x.gap[1]}{x.residue}{x.gap[2]}{x.chain}{x.gap[3]}{x.num}{x.gap[4]}{x.info}\n")
                    x.idx = new_index
                    new_index += 1
                base_index = new_index-2
            elif len(j) == 3 and j[2].atom_info.atom == "N":
                c, o, n = j
                con = [c.atom_info, o.atom_info, n.atom_info]
                file.write("BRANCH" + " " * max(0, (4 - len(str(base_index)))) + f"{base_index}" + " " * max(0, (
                        4 - len(str(base_index)))) + f"{new_index}\n")
                endbranch_list.append((base_index, new_index))
                for x in con:
                    file.write("ATOM" + " " * max(0, (7 - len(
                        str(new_index)))) + f"{new_index}{x.gap[0]}{x.atomnum1}{x.atom}{x.atomtip}{x.atomnum2}{x.gap[1]}{x.residue}{x.gap[2]}{x.chain}{x.gap[3]}{x.num}{x.gap[4]}{x.info}\n")
                    x.idx = new_index
                    new_index += 1
                base_index = new_index - 2
            elif len(j) == 3 and j[2].atom_info.atom == "O":
                c, o, oxt = j
                coo = [c.atom_info, o.atom_info, oxt.atom_info]
                file.write("BRANCH" + " " * max(0, (4 - len(str(base_index)))) + f"{base_index}" + " " * max(0, (
                            4 - len(str(base_index)))) + f"{new_index}\n")
                endbranch_list.append((base_index, new_index))
                for x in coo:
                    file.write("ATOM" + " " * max(0, (7 - len(
                        str(new_index)))) + f"{new_index}{x.gap[0]}{x.atomnum1}{x.atom}{x.atomtip}{x.atomnum2}{x.gap[1]}{x.residue}{x.gap[2]}{x.chain}{x.gap[3]}{x.num}{x.gap[4]}{x.info}\n")
                    x.idx = new_index
                    new_index += 1
                base_index = new_index - 2
            else:
                continue

        else:
            continue
        endbranch_list.sort(reverse=True)
    for j in endbranch_list:
        file.write("ENDBRANCH"+" "*max(0,4-len(str(j[0])))+f"{j[0]}"+" "*max(0, (4 - len(str(base_index))))+f"{j[1]}\n")
    return new_index

def write_last(residue_info_list, file, index):
    new_index = index
    for i in residue_info_list:
        if i["ca"].atom_info.root and i["ca"].atom_info.sidechain:
            new_index = write_sidechain(i["ca"].atom_info.idx, i["ca"], file, new_index)
        else:
            continue
    return new_index

with open(output_pdbqt, "w") as f:
    f.write("ROOT\n")
    f.writelines(root)
    f.write("ENDROOT\n")
    index = write_tail(branch_index, tail, f, index)
    index = write_last(residue_info_list,f, index)

torsdof = 0
with open(output_pdbqt, "r") as f:
    for line in f:
        if line.startswith("BRANCH"):
            torsdof += 1
        else:
            continue

with open(output_pdbqt, "a") as f:
    f.write(f"TORSDOF {torsdof}\n")







