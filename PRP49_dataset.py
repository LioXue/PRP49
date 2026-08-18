import re
import os
import csv

data = r"C:\autodockvina\prp49\outdata"
proteintsv = r"C:\autodockvina\prp49\prodata\pro520.tsv"

def read_atom_info(line):
    pattern = r'^ATOM\s+\d+\s+\d*([A-Za-z]+)\d*\s+(\w+)\s+\w*\s+(\d+)\s+.*$'
    match = re.match(pattern, line)
    if match:
        atom, residue, resnum = match.groups()
        return atom, residue, int(resnum)
    else:
        return None

data_dict = {} #{filename:[MET,ASP]}

def read_pdbqt(path, filename):
    with open(path, "r") as f:
        file_list = [None]*51
        idx = 0
        for line in f:
            if line.startswith("TORSDOF") and isinstance(energy, float):
                file_list[-1] = idx
                data_dict[filename] = file_list
                return energy
            elif line.startswith("ATOM"):
                atom, residue, resnum = read_atom_info(line)
                if atom == "CA":
                    file_list[resnum-1] = residue
                    idx += 1
                else:
                    continue
            elif line.startswith("REMARK VINA RESULT"):
                match = re.search(r'-?\d+\.?\d*', line)
                if match:
                    energy = float(match.group())
                else:
                    energy = None

mapping =  {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
        'ASX': 'D', 'GLX': 'E'
            }

def find_seq(x):
    temp_list = []
    for i in x:
        if not i and len(temp_list) == x[-1]:
            seq = ''.join(temp_list)
            return seq
        elif not i and len(temp_list) != x[-1]:
            return None
        else:
            temp_list.append(mapping[i])

energy_dict = {} #{filename:energy}
lassoseq_dict = {} #{filename:seq}
for filename in os.listdir(data):
    if filename.endswith('.pdbqt'):
        filepath = os.path.join(data, filename)
        try:
            energy_dict[filename] = read_pdbqt(filepath, filename)
            lassoseq_dict[filename] = find_seq(data_dict[filename])
            print(f"已处理: {filename}, energy is {energy_dict[filename]}, lasso seq is {lassoseq_dict[filename]}")
        except Exception as e:
            print(f"处理文件 {filename} 时出错: {e}")
            continue

proteinseq_dict = {} #{filename:seq}
def record_protein(name):
    if name in proteinseq_dict:
        return
    pattern = r'^(\w+)vs(\w+).*$'
    match = re.match(pattern, name)
    if match:
        protein = match.group(1)
        with open(proteintsv, "r", encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if any(protein in str(val) for val in row.values()):
                    proteinseq_dict[name] = row["Sequence"]
                    return
                else:
                    continue
    else:
        print(f'{name}未找到')
        return

for i in data_dict:
    record_protein(i)
    print (f"已处理: {i}, protein seq is {proteinseq_dict[i]}")

keys = energy_dict.keys()

outcsv = r'C:\autodockvina\prp49\outdata\out.csv'
with open(outcsv, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['filename','protein seq', 'lasso seq', 'energy'])
    writer.writeheader()
    for k in keys:
        if proteinseq_dict[k] and lassoseq_dict[k] and energy_dict[k]:
            row = {'filename': k, 'protein seq': proteinseq_dict[k], 'lasso seq': lassoseq_dict[k], 'energy': energy_dict[k]}
            writer.writerow(row)
        else:
            print(f'处理{k}时出错，已跳过')
            continue



