import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import mplcursors 
from matplotlib.lines import Line2D


path_params = "src/distribution_matching/simulation/data/distribution_matching_params.csv"
path_data = "src/distribution_matching/simulation/data/distribution_matching_sim.csv"


ROW_MERGE_COLOR = "#0000d5ff"
CRC_COLOR = "#7600d5ff"
NO_PRECODING_COLOR = "#d00000ff"
PC_COLOR ="#007879ff"

prec_and_colors = {
    "None" : NO_PRECODING_COLOR,
    "RowMerge" : ROW_MERGE_COLOR,
    "ParCheck" : PC_COLOR,
    "CRC" : CRC_COLOR
}

LS_32 = "^"
LS_8 = "<"
LS_1 =  "D"
LS_1_RR = "*"

list_and_marker = {
    32 : LS_32,
    8 : LS_8,
    1 : LS_1,
    0 : LS_1_RR,
    64 : "*"
}

def get_point_for_id(pid):
    df = pd.read_csv(path_data)
    df = df[(df.id == pid)]
    row = df.iloc[0]

    return row["average_distribution"]

def create_dic(s:str) -> dict:
    s = "code=" + s
    s = s.replace(' ', '')
    ls = []
    for f in s.split(')'):
        for part in f.split('('):
            for ppart in f.split(';'):
                ls.append(ppart)
    ls[0] = ls[0].split('(')[1]
    ls[-1] = ls[-1].strip(')')
    return dict(item.split('=') for item in ls if len(item.split('=')) == 2)

def get_legend_for_id(pid):
    df = pd.read_csv(path_params)
    df = df[(df.id == pid)]
    row = df.iloc[0]
    dic = create_dic(row["distribution_matcher"])
    eL = int(dic["encoding_list"])
    fprec = row["distribution_matcher"].split("dynamic_shaping_bit_f=")[1].split("(")[0]
    if fprec == "PCconstraint":
        fprec = "ParCheck"
    rr = int(dic["randomized_rounding"])
    eL -= rr
    return (prec_and_colors[fprec], list_and_marker[eL])

def get_legend():
    labels = ["Precoding"]
    el = [Line2D([0],[0], color="white", alpha=0,lw=3,label="Precoding")]
    for k in prec_and_colors.keys():
        el.append(Line2D([0], [0], color=prec_and_colors[k], lw=3,label=k))
        labels.append(k)
    el.append(Line2D([0],[0],color="white"))
    labels.append("Encoding List Size and RR")
    for k in list_and_marker.keys():
        el.append(Line2D([0],[0],color="black",marker=list_and_marker[k],markerfacecolor="black",markersize=10,lw=0,label=f"{k}"))
        if k != 0:
            labels.append(f"{k}")
        else:
            labels.append(f"{1} with RR")
    for i in range(1):
        el.append(Line2D([0],[0],color="white"))
        labels.append("")
    el.append(Line2D([0], [0], linestyle="--", label="Lower bound", color="darkolivegreen", lw=3))
    labels.append("Rate-Distortion Bound")
    for i in range(3):
        el.append(Line2D([0],[0],color="white"))
        labels.append("")
    return el, labels, 3


def set_legend(fig, fs=20):
    el, lab, cols = get_legend()
    fig.legend(
        handles=el, 
        labels=lab,
        #loc="lower center", 
        #bbox_to_anchor=(0.5, 0.02), 
        ncol=cols, 
        fontsize=fs, 
        title_fontsize=fs, 
        columnspacing=2.5,
        handletextpad=0.5, 
    )

def h2(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))

def get_rate_distortion_bound():
    ds = [i/128 for i in range(0,65)]
    rs = [1 - h2(d) for d in ds]
    return rs, ds

def plot_all_points_pc(N=128,p=0.2396279751674406,fs=20):
    fig = plt.figure()
    df = pd.read_csv(path_params)
    for row in df.iterrows():
        row = row[1]
        if not row["distribution_matcher"].startswith("PolarCode"):
            print(f"unmatching type of code")
            continue
        if row["target_distribution"] != p:
            print(f"unmatching p {p}")
            continue
        dic = create_dic(row["distribution_matcher"])
        n = int(row["blocklength"])
        if n != N:
            print(f"unmatching block length")
            continue
        r = row["info_rate"]
        d = get_point_for_id(row["id"])

        color, mark = get_legend_for_id(row["id"])
        if color == NO_PRECODING_COLOR:
            pass
        plt.scatter(r, d, color=color, marker=mark, s=3*fs, edgecolors="black", linewidths=1)
    rs, ds = get_rate_distortion_bound()
    plt.plot([1 - r for r in rs], ds, linestyle="--", label="Rate-Distortion bound", color="darkolivegreen")

    plt.grid(True, which="both")
    plt.tick_params(axis="both",labelsize=fs)
    plt.ylabel(rf"$D$", fontsize=2*fs)
    plt.xlabel(rf"$R$", fontsize=2*fs)
    set_legend(fig, fs=fs)
    plt.show()

def plot_none_preco(N=128,p=0.2396279751674406,L=1):
    df = pd.read_csv(path_params)
    tuples = []
    print(f"Lower bound for L={L}")
    for row in df.iterrows():
        row = row[1]
        if not row["distribution_matcher"].startswith("PolarCode"):
            #print(f"unmatching type of code")
            continue
        if row["target_distribution"] != p:
            #print(f"unmatching p {p}")
            continue
        dic = create_dic(row["distribution_matcher"])
        n = int(row["blocklength"])
        if n != N:
            #print(f"unmatching block length")
            continue
        r = row["info_rate"]
        d = get_point_for_id(row["id"])
        el = int(dic["encoding_list"])
        if el == 32:
            pass#print(f"found id {row["id"]} with list size {32}")
        color, mark = get_legend_for_id(row["id"])
        if color != NO_PRECODING_COLOR:
            continue
        if mark != list_and_marker[L]:
            continue
        tuples.append((r, d))

        #plt.scatter(r, d, color=color, marker=mark, s=3*fs, edgecolors="black", linewidths=1)
    tup = tuples
    mark = list_and_marker[L]
    tup = sorted(tup,key=lambda t : t[0])
    rs = []
    ps = []
    for t in tup:
        rs.append(t[0])
        ps.append(t[1])
    plt.plot(rs, ps, marker=mark, color=NO_PRECODING_COLOR, label=f"PCDM SCL-{L}")

def get_crc_width(code):
    return int(code.split("width=")[1].split(";")[0])
def get_n_pc_points(code):
    return len(code.split("check_points=[")[1].split("];")[0].replace(" ", "").split(","))

def plot_only_for_pc(N=128,p=0.2396279751674406,fs=20,tcolor=PC_COLOR):
    fig = plt.figure()
    df = pd.read_csv(path_params)
    dico = {}
    count = {}
    for row in df.iterrows():
        row = row[1]
        if not row["distribution_matcher"].startswith("PolarCode"):
            print(f"unmatching type of code")
            continue
        if row["target_distribution"] != p:
            print(f"unmatching p {p}")
            continue
        dic = create_dic(row["distribution_matcher"])
        n = int(row["blocklength"])
        if n != N:
            print(f"unmatching block length")
            continue
        if "PCconstraint" not in row["distribution_matcher"]:
            continue
        r = row["info_rate"]
        d = get_point_for_id(row["id"])
        el = int(dic["encoding_list"])
        rr = int(dic["randomized_rounding"])
        el -= rr
        if el == 0:
            continue
        l = get_n_pc_points(row["distribution_matcher"])
        rv = int(dic["randv"]) if "randv" in dic else 0
        if l != 4:
            continue 
        if rv == 1 and l == 10:
            continue
        k = (el, l, rv)
            
        k2 = (el, l, rv, r)
        color, mark = get_legend_for_id(row["id"])
        if color != tcolor:
            continue
        if k not in dico.keys():
            dico[k] = {r : d}
            count[k2] = 1
        else:
            lis = dico[k]
            if r in lis:
                lis[r] += d
                count[k2] += 1
            else:
                lis[r] = d
                count[k2] = 1
            dico[k] = lis
    for kk in count.keys():
        el, l, rv, r = kk
        k = (el, l, rv)
        lis = dico[k]
        lis[r] /= count[kk]
        dico[k] = lis

        #plt.scatter(r, d, color=color, marker=mark, s=3*fs, edgecolors="black", linewidths=1)
    keys = dico.keys()#[(1,10), (8, 10), (32, 10), (1,4), (8,4), (32,4)]
    for k in keys:
        tup = dico[k]
        el, l, v = k
        mark = list_and_marker[el]
        dd = dict(sorted(tup.items()))
        rs = []
        ps = []
        for t in dd.items():
            rs.append(t[0])
            ps.append(t[1])
        if v == 0:
            color = tcolor
        elif v == 1:
            color = CRC_COLOR
        plt.plot(rs, ps, marker=mark, color=color, label=f"{l}-ParCheck-PCDM SCL-{el}{" (random)" if v == 1 else ""}")
    plot_none_preco(N=N,p=p,L=1)
    rs, ds = get_rate_distortion_bound()
    plt.plot([1 - r for r in rs], ds, linestyle="--", label="Lower bound", color="darkolivegreen")


    plt.grid(True, which="both")
    plt.tick_params(axis="both",labelsize=fs)
    plt.ylabel(rf"$\bar{{p}}$", fontsize=2*fs)
    plt.xlabel(rf"$R$", fontsize=2*fs)
    plt.legend(fontsize=fs)
    plt.show()


def plot_only_for_crc(N=128,p=0.2396279751674406,fs=20,tcolor=CRC_COLOR):
    fig = plt.figure()
    df = pd.read_csv(path_params)
    dico = {}
    for row in df.iterrows():
        row = row[1]
        if not row["distribution_matcher"].startswith("PolarCode"):
            print(f"unmatching type of code")
            continue
        if row["target_distribution"] != p:
            print(f"unmatching p {p}")
            continue
        dic = create_dic(row["distribution_matcher"])
        n = int(row["blocklength"])
        if n != N:
            print(f"unmatching block length")
            continue
        if "CRC" not in row["distribution_matcher"]:
            continue
        r = row["info_rate"]
        d = get_point_for_id(row["id"])
        el = int(dic["encoding_list"])
        w = get_crc_width(row["distribution_matcher"])
        k = (el, w)
            
        color, mark = get_legend_for_id(row["id"])
        if color != tcolor:
            continue
        if k not in dico.keys():
            dico[k] = [(r, d)]
        else:
            dico[k].append((r,d))

        #plt.scatter(r, d, color=color, marker=mark, s=3*fs, edgecolors="black", linewidths=1)
    key = [(1, 11),(1,5), (8,11), (8,5), (32,11), (32, 5), (64,5)]#sorted(dico.keys(), key=lambda t : -t[1])
    for k in key:
        tup = dico[k]
        el, w = k
        mark = list_and_marker[el]
        tup = sorted(tup,key=lambda t : t[0])
        rs = []
        ps = []
        for t in tup:
            rs.append(t[0])
            ps.append(t[1])
        if w == 5:
            color = tcolor
        elif w == 11:
            color = PC_COLOR
        plt.plot(rs, ps, marker=mark, color=color, label=f"{w}-CRC-PCDM SCL-{el}")
    plot_none_preco(N=N,p=p,L=1)
    rs, ds = get_rate_distortion_bound()
    plt.plot([1 - r for r in rs], ds, linestyle="--", label="Lower bound", color="darkolivegreen")


    plt.grid(True, which="both")
    plt.tick_params(axis="both",labelsize=fs)
    plt.ylabel(rf"$\bar{{p}}$", fontsize=2*fs)
    plt.xlabel(rf"$R$", fontsize=2*fs)
    plt.legend(fontsize=fs)
    plt.show()

def get_n_merges(code):
    return int(code.split("n_merges=")[1].split(";")[0])

def plot_only_for_rm(N=128,p=0.2396279751674406,fs=20,tcolor=CRC_COLOR):
    fig = plt.figure()
    df = pd.read_csv(path_params)
    dico = {}
    for row in df.iterrows():
        row = row[1]
        if not row["distribution_matcher"].startswith("PolarCode"):
            #print(f"unmatching type of code")
            continue
        if row["target_distribution"] != p:
            #print(f"unmatching p {p}")
            continue
        dic = create_dic(row["distribution_matcher"])
        n = int(row["blocklength"])
        if n != N:
            #print(f"unmatching block length")
            continue
        if "RowMerge" not in row["distribution_matcher"]:
            #print("no show merge")
            continue
        print(f"found id {row["id"]}")
        r = row["info_rate"]
        d = get_point_for_id(row["id"])
        el = int(dic["encoding_list"])
        rr = int(dic["randomized_rounding"])
        el -= rr
        if el == 0:
            continue
        w = get_n_merges(row["distribution_matcher"])
        k = (el)

        if el == 64:
            continue
            
        color, mark = get_legend_for_id(row["id"])
        if color != ROW_MERGE_COLOR:
            continue
        if k not in dico.keys():
            dico[k] = [(r, d)]
        else:
            dico[k].append((r,d))

        #plt.scatter(r, d, alpha=0.0, label=f"{row["id"]} : n_merges {w}, rate = {r}")

        #plt.scatter(r, d, color=color, marker=mark, s=3*fs, edgecolors="black", linewidths=1)
    key = dico.keys()#[(1, 11),(1,5), (8,11), (8,5), (32,11), (32, 5), (64,5)]#sorted(dico.keys(), key=lambda t : -t[1])
    for k in key:
        tup = dico[k]
        el = k
        mark = list_and_marker[el]
        tup = sorted(tup,key=lambda t : t[0])
        rs = []
        ps = []
        for t in tup:
            rs.append(t[0])
            ps.append(t[1])
        if w == 5:
            color = tcolor
        elif w == 11:
            color = PC_COLOR
        plt.plot(rs, ps, marker=mark, color=color, label=f"RowMerge-PCDM SCL-{el}")
    plot_none_preco(N=N,p=p,L=1)
    plot_none_preco(N=N,p=p,L=32)
    rs, ds = get_rate_distortion_bound()
    plt.plot([1 - r for r in rs], ds, linestyle="--", label="Lower bound", color="darkolivegreen")


    plt.grid(True, which="both")
    #mplcursors.cursor(hover=True)
    plt.tick_params(axis="both",labelsize=fs)
    plt.ylabel(rf"$\bar{{p}}$", fontsize=2*fs)
    plt.xlabel(rf"$R$", fontsize=2*fs)
    plt.legend(fontsize=fs)
    plt.show()



def kldiv(p1a, p1ideal):
    return p1a*np.log2(p1a/p1ideal) + (1-p1a)*np.log2((1-p1a)/(1-p1ideal))

def nom_kl(r, p, ptarg):
    return h2(p) - r +  kldiv(p, ptarg)

def plot_rr(p=0.2396279751674406,fs=20):    
    df = pd.read_csv(path_params)
    for row in df.iterrows():
        row = row[1]
        if not row["distribution_matcher"].startswith("PolarCode"):
            print(f"unmatching type of code")
            continue
        if row["target_distribution"] != p:
            print(f"unmatching p {p}")
            continue
        dic = create_dic(row["distribution_matcher"])
        n = int(row["blocklength"])
        r = row["info_rate"]
        if r*n != int(n*(h2(p) - 1/16)):
            continue
        d = get_point_for_id(row["id"])
        color, mark = get_legend_for_id(row["id"])
        if color != NO_PRECODING_COLOR or mark != LS_1_RR:
            continue
        plt.scatter(n, nom_kl(r, d, p))
    ns = [1 << n for n in range(12)]
    up = [2**(-(n**0.25)) for n in ns]
    plt.plot(ns, up)
    up2 = [np.log2(n)/n for n in ns]
    plt.plot(ns, up2)
    plt.yscale("log")
    plt.xscale("log")
    plt.show()
#plot_rr()
#plot_only_one_type()
#plot_only_for_rm(fs=20)
plot_only_for_pc(fs=20)
#plot_only_for_crc(fs=20)
#plot_all_points_pc(N=512)