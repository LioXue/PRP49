def find_cx(dict):
    p = dict["cg"]+dict["cd"]+dict["ce"]+dict["cz"]
    return p

def a_b(dict):
    """找到dict中ca的原子节点，在sidechain list里添加dict中cb的原子节点位置"""
    dict["ca"].atom_info.sidechain.append(dict["cb"])

def b_g(dict):
    if dict["cb"] is not None:
        dict["cb"].atom_info.sidechain.extend(dict["cg"])
    else:
        pass

def g_d(dict):
    if len(dict["cg"]) != 0:
        dict["cg"][0].atom_info.sidechain.extend(dict["cd"])
    else:
        pass

def d_e(dict):
    if len(dict["cd"]) != 0:
        dict["cd"][0].atom_info.sidechain.extend(dict["ce"])
    else:
        pass

def b_ring(dict):
    """将杂环连接在cb"""
    p = find_cx(dict)
    if dict["cb"] is not None:
        dict["cb"].atom_info.sidechain.extend(p+dict["nx"]+dict["hx"])
    else:
        pass

def b_ring_tyr(dict):
    """Tyr特化的杂环连接"""
    p = find_cx(dict)
    if dict["cb"] is not None:
        dict["cb"].atom_info.sidechain.extend(p)
    else:
        pass

def b_oh(dict):
    """将羟基连接在cb"""
    if dict["cb"] is not None:
        dict["cb"].atom_info.sidechain.extend(dict["ox"]+dict["hx"])
    else:
        pass

def z_oh(dict):
    """将羟基连接在cz"""
    if dict["cz"]:
        dict["cz"][0].atom_info.sidechain.extend(dict["ox"]+dict["hx"])
    else:
        pass

def a_bg(dict):
    """Val，将异丙基连接在ca"""
    dict["ca"].atom_info.sidechain.extend([dict["cb"]]+dict["cg"])

def b_gd(dict):
    """Leu，将异丙基连接在cb"""
    if dict["cb"] is not None:
        dict["cb"].atom_info.sidechain.extend(dict["cg"]+dict["cd"])
    else:
        pass

def b_s(dict):
    """Cys，将巯基连接在cb"""
    if dict["cb"] is not None:
        dict["cb"].atom_info.sidechain.extend(dict["s"]+dict["hx"])
    else:
        pass

def g_sme(dict):
    """Met，将硫-甲基连接在cg"""
    if len(dict["cg"]) != 0:
        dict["cg"][0].atom_info.sidechain.extend(dict["s"]+dict["hx"])
    else:
        pass

def b_cox(dict):
    """Asp/Asn，将羧酸/酰胺连接在cb"""
    if dict["cb"] != None:
        dict["cb"].atom_info.sidechain.extend(dict["cg"]+dict["ox"]+dict["hx"]+dict["nx"])
    else:
        pass

def g_cox(dict):
    """Glu/Gln，将羧酸/酰胺连接在cg"""
    if len(dict["cg"]) != 0:
        dict["cg"][0].atom_info.sidechain.extend(dict["cd"]+dict["ox"]+dict["hx"]+dict["nx"])



def sidechain_find(dict):
    if dict["ca"].atom_info.residue in ("GLY", "PRO"):
        pass
    elif dict["ca"].atom_info.residue == "ALA":
        a_b(dict)
    elif dict["ca"].atom_info.residue == "VAL":
        a_bg(dict)
    elif dict["ca"].atom_info.residue == "LEU":
        a_b(dict)
        b_gd(dict)
    elif dict["ca"].atom_info.residue == "ILE":
        for i in dict["cg"]:
            if i.atom_info.atomnum2 == "2":
                dict["ca"].atom_info.sidechain.extend([dict["cb"],i])
            elif i.atom_info.atomnum2 == "1":
                dict["cb"].atom_info.sidechain.extend(dict["cd"]+[i])
            else:
                pass
    elif dict["ca"].atom_info.residue == "CYS":
        a_b(dict)
        b_s(dict)
    elif dict["ca"].atom_info.residue == "MET":
        a_b(dict)
        b_g(dict)
        g_sme(dict)
    elif dict["ca"].atom_info.residue == "THR":
        dict["ca"].atom_info.sidechain.extend([dict["cb"]]+dict["cg"])
        b_oh(dict)
    elif dict["ca"].atom_info.residue == "SER":
        a_b(dict)
        b_oh(dict)
    elif dict["ca"].atom_info.residue in ("ASN","ASP"):
        a_b(dict)
        b_cox(dict)
    elif dict["ca"].atom_info.residue in ("GLU","GLN"):
        a_b(dict)
        b_g(dict)
        g_cox(dict)
    elif dict["ca"].atom_info.residue == "LYS":
        a_b(dict)
        b_g(dict)
        g_d(dict)
        d_e(dict)
        if len(dict["ce"]) != 0:
            dict["ce"][0].atom_info.sidechain.extend(dict["nx"]+dict["hx"])
        else:
            pass
    elif dict["ca"].atom_info.residue == "ARG":
        a_b(dict)
        b_g(dict)
        g_d(dict)
        if len(dict["cd"]) != 0:
            dict["cd"][0].atom_info.sidechain.extend(dict["nx"]+dict["hx"]+dict["cz"])
        else:
            pass
    elif dict["ca"].atom_info.residue in ("PHE", "TRP", "HIS"):
        a_b(dict)
        b_ring(dict)
    elif dict["ca"].atom_info.residue == "TYR":
        a_b(dict)
        b_ring_tyr(dict)
        z_oh(dict)
    else:
        pass