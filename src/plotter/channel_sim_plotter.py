import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import mplcursors 
from matplotlib.lines import Line2D

path_params = "src/channel_coding/simulation/data/channel_coding_params.csv"
path_sims = "src/channel_coding/simulation/data/channel_coding_sim.csv"
path_dm = "src/channel_coding/simulation/data/distribution_matching_sim.csv"

def intfromtoPlot(toPlot):
    if toPlot == "fer":
        return 5
    elif toPlot == "ber":
        return 6
    else:
        return 0

def plot_polar_ref_ook():
    #ref SCL shaped 3232
    x = [4.9980590062111805, 5.246506211180124, 5.496894409937888, 5.749223602484472, 5.999611801242236, 6.24611801242236]
    y = [0.3485688368472943, 0.13624744228764685, 0.043664920118040094, 0.007091633267391836, 0.001042898691885092, 0.0001306430273120515]
    #plt.plot(x, y, marker="o", label=f"Wiegart : SCL32(enc&dec) shaped")

    #ref SCL unshaped
    x = [4.994177018633541, 5.244565217391305, 5.494953416149068, 5.749223602484472, 5.999611801242236, 6.2519409937888195, 6.498447204968944, 6.748835403726708, 6.997282608695652]
    y = [1, 0.9335750868773576, 0.74241504827435, 0.45538088072371047, 0.20579406028697805, 0.05792318105393938, 0.011561608045966796, 0.0017936237390528277, 0.00021792570404067488]
    plt.plot(x, y, marker="*", label=f"Wiegart : SCL32 unshaped")

def plotRefRM128():
    x_rcub = [
        2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
        3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
        4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0
    ]

    y_rcub = [
        0.0060703,0.0041523,0.0027900,0.0018406,0.0011916,
        0.0007568,0.0004713,0.0002877,0.0001720,0.0001008,
        5.784e-05,3.249e-05,1.787e-05,9.623e-06,5.074e-06,
        2.620e-06,1.326e-06,6.588e-07,3.213e-07,1.541e-07,
        7.289e-08,3.559e-08,1.873e-08,9.264e-09,4.409e-09,
        2.047e-09,9.354e-10,4.229e-10,1.900e-10,8.518e-11,
        3.817e-11
    ]
    x_meta = np.round(np.arange(2.0, 5.0, 0.1), 1)
    x_meta = np.append(x_meta, 4.9)

    y_meta = [
        0.00305257,0.00203604,0.00133077,0.00085181,
        0.00053360,0.00032692,0.00019575,0.00011448,
        6.534e-05,3.637e-05,1.973e-05,1.042e-05,
        5.358e-06,2.678e-06,1.300e-06,6.130e-07,
        2.802e-07,1.241e-07,5.326e-08,2.210e-08,
        8.865e-09,3.433e-09,1.282e-09,4.617e-10,
        1.600e-10,5.335e-11,1.708e-11,5.253e-12,
        1.548e-12,4.369e-13,1.179e-13
    ]

    x_scl4 = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    y_scl4 = [
        0.2295684,0.1090988,0.0401074,0.0111630,
        0.0019682,0.0002692,2.436e-05,1.214e-06,3.969e-08
    ]

    x_scl8 = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    y_scl8 = [
        0.1676133,0.0712103,0.0213935,0.0045071,
        0.0006774,6.725e-05,4.409e-06,2.103e-07,1.099e-08
    ]

    shift = 10 * np.log10(2*0.46875)
    plt.semilogy(x_rcub + shift, y_rcub, label="Random coding union bound")
    plt.semilogy(x_meta + shift, y_meta, label="Meta converse bound")

    plt.plot(x_scl4 + shift, y_scl4, marker="o", linestyle="-", label="Zunker : RowMerged SCL-4")
    plt.plot(x_scl8 + shift, y_scl8, marker="o", linestyle="-", label="Zunker : RowMerged SCL-8")

    plt.xlim(2,5)
    plt.ylim(1/(10**8), 1)

def plotRefRM256():
    x_rcub = np.arange(1.0, 4.01, 0.1)

    y_rcub = [
        0.0097707,0.0065552,0.0043018,0.0027591,0.0017283,
        0.0010565,0.0006298,0.0003658,0.0002068,0.0001137,
        6.084e-05,3.159e-05,1.592e-05,7.783e-06,3.685e-06,
        1.690e-06,7.497e-07,3.215e-07,1.332e-07,5.329e-08,
        2.057e-08,7.659e-09,2.749e-09,9.510e-10,3.171e-10,
        1.019e-10,3.159e-11,9.452e-12,2.733e-12,7.651e-13,
        2.079e-13
    ]

    x_meta = np.arange(1.0, 3.51, 0.05)
    x_meta = np.append(x_meta, [3.6,3.7,3.8,3.9,4.0])

    y_meta = [
        0.0061419,0.0050126,0.0040673,0.0032810,0.0026308,
        0.0020968,0.0016608,0.0013072,0.0010223,0.0007944,
        0.0006132,0.0004701,0.0003580,0.0002707,0.0002033,
        0.0001516,0.0001122,8.246e-05,6.013e-05,4.352e-05,
        3.125e-05,2.227e-05,1.574e-05,1.103e-05,7.676e-06,
        5.294e-06,3.620e-06,2.454e-06,1.649e-06,1.098e-06,
        7.254e-07,4.745e-07,3.076e-07,1.975e-07,1.256e-07,
        7.912e-08,4.934e-08,3.047e-08,1.862e-08,1.127e-08,
        6.750e-09,3.999e-09,2.344e-09,1.359e-09,7.798e-10,
        4.423e-10,2.480e-10,1.375e-10,7.535e-11,4.080e-11,
        2.183e-11
    ]

    x_scl4 = [1.0,1.5,2.0,2.5,3.0,3.5,4.0]
    y_scl4 = [
        0.1177478,0.0387173,0.0092598,
        0.0015507,0.0001640,1.196e-05,5.451e-07
    ]

    x_scl8 = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]
    y_scl8 = [
        0.3615328,0.1932740,0.0721240,
        0.0194867,0.0035448,0.0004158,
        2.848e-05,1.389e-06,3.878e-08
    ]

    plt.semilogy(x_rcub, y_rcub, label="Random coding union bound")
    #plt.semilogy(x_meta, y_meta, label="Meta converse bound")

    plt.semilogy(x_scl4, y_scl4, marker="*", linestyle="-", label="Zunker : RowMerged SCL-4")
    plt.semilogy(x_scl8, y_scl8, marker="o", linestyle="-", label="Zunker : RowMerged SCL-8")

    plt.xlim(2,4)
    plt.ylim(1/(10**8), 1)


def parse_list(lst:str) -> list:
    lst = lst.strip(']')
    lst = lst.strip['[']
    lst = lst.replace(' ', '')
    elems = lst.split(',')
    return elems

def crc_poly_to_str(poly:str) -> str:
    pol = [int(s) for s in parse_list(lst)]
    deg = len(pol) - 1
    s = rf""
    for p in pol:
        if p == 1:
            s += rf"x^{deg}"
        deg -= 1
    return s


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

def find_polar_legen(rate, shaping_rate, s:str) -> str:
    di = create_dic(s)
    encL = int(di["encoding_list"])
    decL = int(di["decoding_list"])
    desc = f"scl{encL}&{decL}"
    froz_prec = di["dynamic_frozen_bit_f"]
    shap_pref = di["dynamic_shaping_bit_f"]




def find_id(blocklength, rate, shaping_rate, constellation, code):
    df = pd.read_csv(path_params)
    df = df[
        (df["blocklength"] == blocklength)
        &(df["rate"] == rate)
        &(df["shaping_rate"] == rate)
        &(df["constellation"] == constellation)
        &(df["code"] == code)
    ]

def plot_for_id(pid, toplot="fer", linestyle="-", mark="o", leg="", color="red", markersize=10, markerfacecolor=None, allowed_buggy=False, show=False):
    dfdm = pd.read_csv(path_dm)
    dfsim = pd.read_csv(path_sims)
    dfdm = dfdm[(dfdm["id"] == pid)]
    if not allowed_buggy:
        dfparam = pd.read_csv(path_params)
        row = dfparam.iloc[pid-1]
        if row["buggy"] == 1:
            return
    if markerfacecolor == None:
        markerfacecolor = color
    dfsim = dfsim[
        (dfsim["id"] == pid)
        &(dfsim["fec"] >= 50)
    ]
    merged = dfsim.merge(dfdm, left_on="snr", right_on="target_snr", how="left")

    x = merged["observed_snr"].fillna(merged["snr"])
    y = merged[toplot]

    plt.plot(x, y, marker=mark, label=leg, linestyle=linestyle, markersize=markersize, color=color, markerfacecolor=markerfacecolor, markeredgecolor="black")
    if show :
        plt.yscale("log")
        plt.show()


def plot_for_ids(pids):
    #AI generated script
    df = pd.read_csv(path_params)
    N = -1
    R = -1
    fs = 20
    fig = plt.figure()
    for pid in pids:
        row = df.iloc[pid - 1]
        if N == -1:
            N = row["blocklength"]
        if N != row["blocklength"]:
            print(f"{pid} is not valid (blocklength != {N})")
            return
        r = row["rate"]
        if R == -1:
            R = r
        if r != R:
            print(f"{pid} is not valid (rate != {R})")
            return

        sr = row["shaping_rate"]
        dic = create_dic(row["code"])
        eL = int(dic["encoding_list"])
        dL = int(dic["decoding_list"])

        fprec = (
            row["code"]
            .split("dynamic_frozen_bit_f=")[1]
            .split("(")[0]
        )
        if fprec == "PCconstraint":
            fprec = "ParCheck"

        sprec = (
            row["code"]
            .split("dynamic_shaping_bit_f=")[1]
            .split("(")[0]
        )
        if sprec == "PCconstraint":
            sprec = "ParCheck"

        if sprec != "ParCheck":
            continue#pass#continue
        
        
        if sprec == "None":
            color = "#d00000ff"
        elif sprec == "RowMerge":
            color = "#0000d5ff"
        elif sprec == "ParCheck":
            color = "#007879ff"
        elif sprec == "CRC":
            color = "#7600d5ff"
        else:
            print(f"sprec = {sprec}")
            color = "black"

        if fprec == "CRC":
            mark = "^"
        elif fprec == "RowMerge":
            mark = "<"
        elif fprec == "ParCheck":
            mark = "D"
        elif fprec == "None":
            mark = "h"
        else:
            print(f"fprec = {fprec}")
            mark = "o"

        if eL == 32:
            lines = "-"
        elif eL == 8:
            lines = "--"
        elif eL == 1:
            lines = ":"
        else:
            print(f"eL = {eL}")
            lines = "-"

       
        if dL == 32:
            markerfacecolor = color
        elif dL == 8:
            markerfacecolor = f"{color[:-2]}88"#f"light{color}"#"white"
        elif dL == 1:
            markerfacecolor = "none"
        else:
            print(f"dL = {dL}")
            markerfacecolor = color

        leg = (
            rf"$R_S$={sr:.2f} "
            rf"{fprec}-precoding "
            rf"{sprec}-shaping "
            rf"({pid})"
        )

        plot_for_id(
            pid,
            leg=leg,
            color=color,
            markersize=int(fs*5/12),
            mark=mark,
            linestyle=lines,
            markerfacecolor=markerfacecolor
        )

   
    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2 * fs)
    plt.xlabel("SNR [dB]", fontsize=2 * fs)
    plt.tick_params(
        axis="both",
        labelsize=fs
    )
    mplcursors.cursor(hover=True)
    # plt.xlim(2, 7)
    #mplcursors.cursor(hover=True)
    plt.yscale("log")
    legend_sprec = [
        Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Shaping"
        ),
        Line2D(
            [0], [0],
            color="#d00000ff",
            lw=3,
            label="None"
        ),
        Line2D(
            [0], [0],
            color="#0000d5ff",
            lw=3,
            label="RowMerge"
        ),
        Line2D(
            [0], [0],
            color="#007879ff",
            lw=3,
            label="ParCheck"
        ),
        Line2D(
            [0], [0],
            color="#7600d5ff",
            lw=3,
            label="CRC"
        ),
    ]
    legend_fprec = [
        Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Precoding"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="^",
            lw=0,
            markersize=10,
            label="CRC"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="<",
            lw=0,
            markersize=10,
            label="RowMerge"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="D",
            lw=0,
            markersize=10,
            label="ParCheck"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="h",
            lw=0,
            markersize=10,
            label="None"
        ),
    ]
    legend_eL = [
         Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Encoding List Size"
        ),
        Line2D(
            [0], [0],
            color="black",
            lw=3,
            linestyle=":",
            label="1"
        ),
        Line2D(
            [0], [0],
            color="black",
            lw=3,
            linestyle="--",
            label="8"
        ),
        Line2D(
            [0], [0],
            color="black",
            lw=3,
            linestyle="-",
            label="32"
        ),
    ]
    legend_dL = [
        Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Deocoding List Size"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="D",
            markerfacecolor="black",
            markersize=10,
            lw=0,
            label="32"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="D",
            markerfacecolor="#00000088",
            markersize=10,
            lw=0,
            label="8"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="D",
            markerfacecolor="none",
            markeredgecolor="black",
            markersize=10,
            lw=0,
            label="1"
        ),
    ]


    legend_elements = [legend_sprec[0], legend_sprec[1], legend_sprec[2], legend_sprec[3], legend_sprec[4], legend_fprec[0], legend_fprec[1], legend_fprec[2], legend_fprec[3], legend_fprec[4], legend_eL[0], legend_eL[1], legend_eL[2], legend_eL[3], None, legend_dL[0], legend_dL[1], legend_dL[2], legend_dL[3], None, ]
    #legend_elements = [legend_sprec[0], legend_fprec[0], legend_eL[0], legend_dL[0], legend_sprec[1], legend_fprec[1], legend_eL[1], legend_dL[1], legend_sprec[2], legend_fprec[2], legend_eL[2], legend_dL[2], legend_sprec[3], legend_fprec[3], None, None, ]
   
    fig.legend(
        handles=legend_elements, 
        labels=["Shaping", "None", "RowMerge", "ParCheck", "CRC", "Precoding", "CRC", "RowMerge", "ParCheck", "None",  "Encoding List size", "1", "8", "32", "", "Decoding List size", "32", "8", "1", "", ],
        #loc="lower center", 
        #bbox_to_anchor=(0.5, 0.02), 
        ncol=4, 
        fontsize=fs, 
        title_fontsize=fs, 
        columnspacing=2.5,
         handletextpad=0.5, )
    # On réserve de la place sous le graphe
    #fig.subplots_adjust(bottom=0.27)
    # sprec
    '''
    fig.legend(
        handles=legend_sprec,
        title="sprec",
        loc="lower center",
        bbox_to_anchor=(0.17, 0.02),
        ncol=1,
        fontsize=fs // 2,
        title_fontsize=fs // 2,
    )
    # fprec
    fig.legend(
        handles=legend_fprec,
        title="fprec",
        loc="lower center",
        bbox_to_anchor=(0.43, 0.02),
        ncol=1,
        fontsize=fs // 2,
        title_fontsize=fs // 2,
    )
    # eL
    fig.legend(
        handles=legend_eL,
        title="eL",
        loc="lower center",
        bbox_to_anchor=(0.67, 0.02),
        ncol=1,
        fontsize=fs // 2,
        title_fontsize=fs // 2,
    )
    # dL
    fig.legend(
        handles=legend_dL,
        title="dL",
        loc="lower center",
        bbox_to_anchor=(0.90, 0.02),
        ncol=1,
        fontsize=fs // 2,
        title_fontsize=fs // 2,
    )
    '''
    plt.show()



def plot_for_ids2(pids):
    df = pd.read_csv(path_params)
    N = -1
    R = -1
    fs = 20
    for pid in pids:
        row = df.iloc[pid-1]
        if N == -1:
            N = row["blocklength"]
        if N != row["blocklength"]:
            print(f"{pid} is not valid (blocklength != {N})")
            return
        N = row["blocklength"]
        r = row["rate"]
        if R == -1:
            R = r
        if r != R:
            print(f"{pid} is not valid (rate != {R})")
            return
        sr = row["shaping_rate"]
        dic = create_dic(row["code"])
        eL = int(dic["encoding_list"])
        dL = int(dic["decoding_list"])
        fprec = row["code"].split("dynamic_frozen_bit_f=")[1].split("(")[0]
        if fprec == "PCconstraint":
            fprec = "ParCheck"
        sprec = row["code"].split("dynamic_shaping_bit_f=")[1].split("(")[0]
        if sprec == "PCconstraint":
            sprec = "ParCheck"

        if sprec == "None":
            color="red"
        elif sprec == "RowMerge":
            color="blue"
        elif sprec == "ParCheck":
            color="green"
        elif sprec == "CRC":
            color="purple"
        else:
            print(f"sprec = {sprec}")

        if fprec != "CRC":
            pass#continue

        if fprec == "CRC":
            mark = "^"
        elif fprec == "RowMerge":
            mark = "<"
        elif fprec == "ParCheck":
            mark = "D"
        elif fprec == "None":
            mark = "h"
        else:
            print(f"fprec = {fprec}")

        if eL == 32:
            lines = "-"
        elif eL == 8:
            lines = "--"
        elif eL == 1:
            lines = ":"
        else:
            print(f"eL = {eL}")
        if dL == 32:
            mark2 = color
        elif dL == 8:
            mark2 = "white"
        elif dL == 1:
            mark2 = "none"
        else:
            print(f"dL = {dL}")
        leg = rf"$R_S$={sr:.2f} {fprec}-precoding {sprec}-shaping ({pid})"
        plot_for_id(pid, leg=leg, color=color, markersize=fs//2, mark=mark, linestyle=lines, markerfacecolor=mark2)
    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2*fs)
    plt.xlabel("SNR [dB]", fontsize=2*fs)
    plt.tick_params(axis="both",labelsize=fs)
    #plt.xlim(2,7)
    mplcursors.cursor(hover=True)
    plt.yscale("log")
    #plt.legend(fontsize=fs)
    plt.show()

def plot_for_shaping_pc():
    df = pd.read_csv(path_params)
    pids = range(154,189)
    N = 64
    R = 0.5
    fs = 10
    for kk in [4,6,8,10,16]:
        for pid in pids:
            row = df.iloc[pid-1]
            if N != row["blocklength"]:
                print(f"{pid} is not valid (blocklength != {N})")
                return
            N = row["blocklength"]
            r = row["rate"]
            if r != R:
                print(f"{pid} is not valid (rate != {R})")
                return
            sr = row["shaping_rate"]
            if not 0.77 <= sr <= 0.85:
                continue
            dic = create_dic(row["code"])
            eL = dic["encoding_list"]
            dL = dic["decoding_list"]
            k = len(row["code"].split("dynamic_shaping_bit_f=PCconstraint(check_points=")[1].split(";")[0].replace(" ", "").strip("[").strip("]").split(","))
            #if k not in [kk]:
            #    continue
            fprec = row["code"].split("dynamic_frozen_bit_f=")[1].split("(")[0]
            leg = rf"$R_S$={sr:.2f} SCL{eL}&{dL} {fprec}-precoding {k}-pc-check"
            plot_for_id(pid, leg=leg)
        plt.grid(True, which="both")
        plt.ylabel("FER", fontsize=fs)
        plt.xlabel("SNR [dB]", fontsize=fs)
        plt.xlim(2,7)
        plt.yscale("log")
        plt.legend(fontsize=fs)
        plt.show()

def plot_5g_comp_L32_CRC11(toPlot="fer"):
    #ids = [48, 49, 50, 51]
    #ids = [52, 53, 54]
    #ids = [55]
    ids = [58, 59]
    df = pd.read_csv(path_params)
    for pid in ids:
        plot_for_id(pid, toPlot, leg="CRC-11 SCL32")
        row = df.iloc[pid-1]
        data = np.loadtxt(f"data/polar_sc/results_polar_R={0.5}_M={2}_n={row["blocklength"]}_L={32}_CRC-{11}.txt", skiprows=1)
        plt.plot(data[:,0], data[:,intfromtoPlot(toPlot)] / data[:,4], marker='*', label=f"5G code N={row["blocklength"]}, R={0.5},L={32},CRC={11}")
        plt.xlabel("SNR [dB]")
        plt.yscale("log")
        plt.grid(True, which="both")
        plt.title(f"5G with SCL32 and CRC11, comparison")
        plt.legend()
        plt.show()

def get_avg_shaping_perf(pid):
    df = pd.read_csv(path_dm)
    df = df[(df["id"] == pid)]
    p1 = sum(df["observed_dist"]*df["frames"]) / sum(df["frames"])
    return p1

def kldiv(p1a, p1ideal):
    return p1a*np.log2(p1a/p1ideal) + (1-p1a)*np.log2((1-p1a)/(1-p1ideal))


def h2(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))

def kldiv(p1a, p1ideal):
    return p1a*np.log2(p1a/p1ideal) + (1-p1a)*np.log2((1-p1a)/(1-p1ideal))

def get_rate_distortion_bound():
    ds = [i/128 for i in range(0,65)]
    rs = [1 - h2(d) for d in ds]
    return rs, ds

def plot_rate_loss(ids):
    df = pd.read_csv(path_params)
    bl = [(pid, df.iloc[pid-1]) for pid in ids]
    d = {}
    for pid, row in bl:
        p1 = get_avg_shaping_perf(pid)
        shap_rate = row["shaping_rate"]
        prec = "CRC" if "dynamic_shaping_bit_f=CRC" in row["code"] else "none"
        dic = create_dic(row["code"])
        encl = int(dic["encoding_list"])
        rl = h2(p1) - shap_rate
        if (shap_rate, prec, encl) not in d:
            d[(shap_rate, prec, encl)] = ([row["blocklength"]], [rl])
        else:
            bls, rls = d[(shap_rate, prec, encl)]
            bls.append(row["blocklength"])
            rls.append(rl)
            d[(shap_rate, prec, encl)] = (bls, rls)
    print(d)
    for k in d.keys():
        bls, rls = d[k]
        plt.scatter(bls, rls, label=f"{k}")
    plt.yscale("log")
    plt.legend()
    plt.show()

def plot_low_shaping():
    #ids = [1,2,3,4,5,6,7,8,9,10,11,12,13] #with reaveal and poly 7
    #ids = [14,15,16,17,18,19,20,21] #with reaveal and poly 3
    #ids = [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33] #without reaveal and poly 3
    #ids = [(1, 14, 23), (2,15,24),(3,16,25), (4,17,26),(5,18,27),(6,19,29),(7,20,30),(8,21,31),(9,32),(10, 33),(11),(12),(13)]
    #ids = [35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47] #without reaveal, poly rev(3)
    ids = [60, 61, 62, 63, 64]
    ids = range(60, 88)
    df = pd.read_csv(path_params)
    for pid in ids:
        #for pid in pids:
        row = df.iloc[pid-1]
        cc = row["channel_changed"]
        fn = "known_at_decoding=[1" in row["code"]
        leg = f"{pid} : CRC-3 SCL32&32 {"(with code unaware of channel)" if cc == 0 else "(with code aware of channel)"} {"crc known at decoding" if fn else "crc unknown at decoding"}"
        plot_for_id(pid, "fer", leg=leg)
        d = int(32*(1 - row["shaping_rate"]))
        data = np.loadtxt(f"data/polar_sc/fer_ook_short_scl_n=5_d={d}.txt", skiprows=1)
        plt.plot(data[:,0], data[:,2],marker="*", label=f"ref curve")
        plt.xlabel("SNR [dB]")
        plt.yscale("log")
        plt.grid(True, which="both")
        plt.title(f"Shaping with {d} bits, comparison")
        plt.legend()
        plt.show()

def plot_shaping_perf3():
    df = pd.read_csv("src/plotter/shaping_sim3.csv")
    p02 = 0.23962797516744064#0.23962797516744064
    p03 = 0.3521832003658291
    fs = 20
    ms = 10
    df021 = df[(df["target_p"] != p03) & (df["L"] == 1)]
    df021 = df021.sort_values("1_shaping_rate") 
    #plt.plot(1 - df021["1_shaping_rate"], df021["obtained_p"], marker="*", label="p=0.24 L=1", color="red", markersize=ms)
    df021 = df[(df["target_p"] != p03) & (df["L"] == 32)]
    df021 = df021.sort_values("1_shaping_rate") 
    plt.plot(1 - df021["1_shaping_rate"], df021["obtained_p"], marker="^", label="p=0.24 L=32", color="#d00000ff", markersize=6)

    rs, ds = get_rate_distortion_bound()
    plt.plot([1 - r for r in rs], ds, linestyle="--", label="Lower bound", color="darkolivegreen")
    
def plot_shaping_perf_for_ids(pids):
    df = pd.read_csv(path_params)
    N = -1
    R = -1
    fs = 20
    fig = plt.figure()
    dico = {}
    for pid in pids:
        row = df.iloc[pid-1]
        if N == -1:
            N = row["blocklength"]
        if N != row["blocklength"]:
            print(f"{pid} is not valid (blocklength != {N})")
            return
        N = row["blocklength"]
        r = row["rate"]
        if R == -1:
            R = r
        if r != R:
            print(f"{pid} is not valid (rate != {R})")
            return
        sr = row["shaping_rate"]
        dic = create_dic(row["code"])
        eL = int(dic["encoding_list"])
        dL = int(dic["decoding_list"])
        fprec = row["code"].split("dynamic_frozen_bit_f=")[1].split("(")[0]
        if fprec == "PCconstraint":
            fprec = "ParCheck"
        sprec = row["code"].split("dynamic_shaping_bit_f=")[1].split("(")[0]
        if sprec == "PCconstraint":
            sprec = "ParCheck"

        key = (sr, eL, sprec)
        if key in dico:
            dico[key] += get_avg_shaping_perf(pid)
            dico[key] /= 2
        else:
            dico[key] = get_avg_shaping_perf(pid)
    for k in dico.keys():
        sr, eL, sprec = k
        if sprec == "None":
            color = "#d00000ff"
        elif sprec == "RowMerge":
            color = "#0000d5ff"
        elif sprec == "ParCheck":
            color = "#007879ff"
        elif sprec == "CRC":
            color = "#7600d5ff"
        else:
            print(f"sprec = {sprec}")
            color = "black"

        if eL == 32:
            mark = "^"
        elif eL == 8:
            mark = "<"
        elif eL == 1:
            mark = "D"

        val = dico[k]

        plt.scatter(k[0], val, marker=mark, color=color, s=3*fs)
    plot_shaping_perf3()
    plt.grid(True, which="both")
    plt.ylim(0,0.55)
    #plt.yscale("log")
    plt.tick_params(axis="both",labelsize=fs)
    plt.ylabel(rf"$\bar{{p}}$", fontsize=2*fs)
    plt.xlabel(rf"$R_M$", fontsize=2*fs)

    legend_elements = [
        Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Shaping"
        ),
        Line2D(
            [0], [0],
            color="#d00000ff",
            lw=3,
            label="None"
        ),
        Line2D(
            [0], [0],
            color="#0000d5ff",
            lw=3,
            label="RowMerge"
        ),
        Line2D(
            [0], [0],
            color="#007879ff",
            lw=3,
            label="ParCheck"
        ),
        Line2D(
            [0], [0],
            color="#7600d5ff",
            lw=3,
            label="CRC"
        ),
         Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Precoding"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="D",
            lw=0,
            markersize=10,
            label="1"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="<",
            lw=0,
            markersize=10,
            label="8"
        ),
        Line2D(
            [0], [0],
            color="black",
            marker="^",
            lw=0,
            markersize=10,
            label="32"
        ),
        Line2D([0], [0], color="#00000000"),
        Line2D(
            [0], [0],
            color="#00000000",
            lw=3,
            label="Lower bounds"
        ),
        Line2D(
            [0], [0],
            color="#d00000ff",
            marker="^",
            markersize=6,
            lw=3,
            label=""
        ),
        Line2D(
            [0], [0],
        linestyle="--", label="Lower bound", color="darkolivegreen", lw=3),
        Line2D([0], [0], color="#00000000"),
        Line2D([0], [0], color="#00000000"),
    ]
    fig.legend(
        handles=legend_elements, 
        labels=["Shaping", "None", "RowMerge", "ParCheck", "CRC", "Encoding List size", "1", "8", "32", "", "", "L32 None-Shaping", "lower bound", "", "",],
        #loc="lower center", 
        #bbox_to_anchor=(0.5, 0.02), 
        ncol=3, 
        fontsize=fs, 
        title_fontsize=fs, 
        columnspacing=2.5,
         handletextpad=0.5, )
    #plt.legend(fontsize=fs)
    plt.show()

def plot_shaping_perf_for_ids2(pids):
    df = pd.read_csv(path_params)
    N = -1
    fs = 20
    ms = 10
    for pid in pids:
        row = df.iloc[pid-1]
        if N == -1:
            N = row["blocklength"]
        if N != row["blocklength"]:
            print(f"{pid} is not valid (blocklength != {N})")
            return
        N = row["blocklength"]
        sr = row["shaping_rate"]
        dic = create_dic(row["code"])
        eL = dic["encoding_list"]
        dL = dic["decoding_list"]
        fprec = row["code"].split("dynamic_frozen_bit_f=")[1].split("(")[0]
        if fprec == "PCconstraint":
            fprec = "ParCheck"
        sprec = row["code"].split("dynamic_shaping_bit_f=")[1].split("(")[0]
        if sprec == "PCconstraint":
            sprec = "ParCheck"
        leg = rf"{fprec}-precoding {sprec}-shaping"
        if sprec == "None":
            color="red"
        elif sprec == "RowMerge":
            color="blue"
        elif sprec == "ParCheck":
            color="green"
        elif sprec == "CRC":
            color="purple"
        plt.scatter(sr, get_avg_shaping_perf(pid), label=leg, marker="+", s=20*fs, color=color)

    plot_shaping_perf3()
    plt.grid(True, which="both")
    plt.ylim(0,0.5)
    #plt.yscale("log")
    plt.tick_params(axis="both",labelsize=fs)
    plt.ylabel(rf"$\bar{{p}}$", fontsize=fs)
    plt.xlabel(rf"$R_S$", fontsize=fs)

    plt.legend(fontsize=fs)
    plt.show()

def get_lower_bound_on_l1_dist_sm3():
    def lower_bound(p_ideal, p_obtained):
        h2_idea = h2(p_ideal)
        h2_obt = h2(p_obtained)
        d = kldiv(p_obtained, p_ideal)
        fact = -1/np.log(p_ideal)
        return fact * (h2_obt - h2_idea + d)

    def upper_bound1(p_ideal, p_obtained, rate, N):
        expon = N * (h2(p_obtained) - rate + kldiv(p_obtained, p_ideal))
        return np.sqrt(1 + np.exp(-expon))

    def upper_bound2(p_ideal, p_obtained, rate, N):
        r = h2(p_obtained) - rate + kldiv(p_obtained, p_ideal)
        return np.sqrt(0.5*N*r)

    df = pd.read_csv("src/plotter/shaping_sim3.csv")
    p02 = 0.23962797516744064#0.23962797516744064
    fs = 20
    ms = 10
    df021 = df[(df["L"] == 1)]
    df021 = df021.sort_values("1_shaping_rate") 
    plt.plot(1 - df021["1_shaping_rate"], [lower_bound(p02, p) for p in df021["obtained_p"]], marker="*", label="p=0.24 L=1", color="red", markersize=ms)
    plt.plot(1 - df021["1_shaping_rate"], [upper_bound2(p02, p, 1-r, 128) for r, p in zip(df021["1_shaping_rate"], df021["obtained_p"])], marker="*", label="p=0.24 L=1", color="red", markersize=ms)
    plt.plot(1 - df021["1_shaping_rate"], [upper_bound1(p02, p, 1-r, 128) for r, p in zip(df021["1_shaping_rate"], df021["obtained_p"])], marker="*", label="p=0.24 L=1", color="red", markersize=ms)

    df021 = df[(df["L"] == 32)]
    df021 = df021.sort_values("1_shaping_rate")
    plt.plot(1 - df021["1_shaping_rate"], [lower_bound(p02, p) for p in df021["obtained_p"]], marker="*", label="p=0.24 L=32", color="green", markersize=ms)
    plt.plot(1 - df021["1_shaping_rate"], [upper_bound2(p02, p, 1-r, 128) for r, p in zip(df021["1_shaping_rate"], df021["obtained_p"])], marker="*", label="p=0.24 L=32", color="green", markersize=ms)
    plt.plot(1 - df021["1_shaping_rate"], [upper_bound1(p02, p, 1-r, 128) for r, p in zip(df021["1_shaping_rate"], df021["obtained_p"])], marker="*", label="p=0.24 L=32", color="green", markersize=ms)

    h22 = h2(p02)
    plt.axvline(h22, linestyle="--", label="h2(0.24)", color="darkred")

    plt.yscale("log")

    plt.legend()
    plt.show()


def plot_smt():
    ids = range(88, 104)
    df = pd.read_csv(path_params)
    for pid in ids:
        row = df.iloc[pid-1]
        rate = row["rate"]
        shap_rate = row["shaping_rate"]
        #if row["shaping_rate"] != 0.765625:
        #    continue
        prec = "CRC11" if "dynamic_frozen_bit_f=CRC" in row["code"] else ""
        shap_prec = "CRC11" if "dynamic_shaping_bit_f=CRC" in row["code"] else ""
        d = create_dic(row["code"])
        c = f"SCL{d["encoding_list"]}&{d["decoding_list"]}"
        if d["encoding_list"] != d["decoding_list"] or int(d["decoding_list"]) not in [4, 8]:
            continue
        leg = f"{shap_prec}&{prec} {c} with rate {rate} and shaping rate {shap_rate}"
        plot_for_id(pid, "fer", leg=leg)
    plt.xlabel("SNR [dB]")
    plt.yscale("log")
    plt.grid(True, which="both")
    plt.title(f"CRC-11 precoding for shaping, blocklength 128, encoding and decoding list size 8, and rate 0.5. Different shaping rate")
    plt.legend()
    plt.show()

def plot_rate_loss2():
    df = pd.read_csv("src/plotter/shaping_sim.csv")
    dic = {}
    p = 0.23962797516744064
    res = 1/16
    ds = [1-h2(p) + i*res for i in range(-2,3)]
    def find_d(r, N):
        for d in ds:
            if d == ds[2] and int(N*d) / N == r:
                return d
        return -1
    for row in df.iterrows():
        row = row[1]
        L = row["L"]
        N = row["N"]
        r = row["1_shaping_rate"]
        p1 = row["obtained_p"]
        key = (L, find_d(r, N))
        if key[1] == -1:
            continue
        if key not in dic:
            dic[key] = {N:h2(p1) - (1-r)}
        else:
            dic[key][N] = h2(p1) - (1-r)
    for k in dic.keys():
        L, r = k
        val = dic[k]
        plt.plot(val.keys(), val.values(), label=f"L={L}, rate={1-r}")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Blocklength")
    plt.ylabel("Rate loss")
    plt.legend()
    plt.show()

def plot_shaping_perf():
    df = pd.read_csv("src/plotter/shaping_sim2.csv")
    p02 = 0.23962797516744064#0.23962797516744064
    p03 = 0.3521832003658291
    fs = 20
    ms = 10
    df021 = df[(df["target_p"] != p03) & (df["L"] == 1)]
    df021 = df021.sort_values("shaping_rate") 
    plt.plot(1 - df021["shaping_rate"], df021["obtained_p"], marker="*", label="p=0.24 L=1", color="red", markersize=ms)
    df021 = df[(df["target_p"] != p03) & (df["L"] == 32)]
    df021 = df021.sort_values("shaping_rate") 
    plt.plot(1 - df021["shaping_rate"], df021["obtained_p"], marker="o", label="p=0.24 L=32", color="red", markersize=ms)
    df021 = df[(df["target_p"] == p03) & (df["L"] == 1)]
    df021 = df021.sort_values("shaping_rate") 
    plt.plot(1 - df021["shaping_rate"], df021["obtained_p"], marker="*", label="p=0.35 L=1", color="blue", markersize=ms)
    df021 = df[(df["target_p"] == p03) & (df["L"] == 32)]
    df021 = df021.sort_values("shaping_rate") 
    plt.plot(1 - df021["shaping_rate"], df021["obtained_p"], marker="o", label="p=0.35 L=32", color="blue", markersize=ms)

    h22 = h2(p02)
    h23 = h2(p03)
    #plt.axvline(h22, linestyle="--", label="h2(0.24)", color="darkred")
    #plt.axvline(h23, linestyle="--", label="h2(0.35)", color="darkblue")
    rs, ds = get_rate_distortion_bound()
    plt.plot([1 - r for r in rs], ds, linestyle="--", label="Lower bound", color="darkolivegreen")
    plt.grid(True, which="both")
    plt.xlim(0.45,0.95)
    plt.ylim(0,0.5)
    #plt.yscale("log")
    plt.tick_params(axis="both",labelsize=30)
    plt.ylabel(rf"$\bar{{p}}$", fontsize=30)
    plt.xlabel(rf"$R$", fontsize=30)

    plt.legend(fontsize=25)
    plt.show()

def plot_normalized_info_div():
    N = 512
    fs = 20
    #D/N = H(~X) - R_s + D(~X, X)
    def nom_kl(r, p, ptarg):
        return h2(p) - r +  kldiv(p, ptarg)
    df = pd.read_csv("src/plotter/shaping_sim2.csv")
    p02 = 0.230962797516744064#0.23962797516744064
    p03 = 0.3521832003658291
    df021 = df[(df["target_p"] != p03) & (df["L"] == 1)]
    df021 = df021.sort_values("shaping_rate")
    shaping_rates = 1 -df021["shaping_rate"]
    res = [nom_kl(r, p, p02) for r,p in zip(shaping_rates, df021["obtained_p"])]
    plt.plot(shaping_rates, res, label="p=0.24 L=1", marker="*", color="orangered")
    df021 = df[(df["target_p"] != p03) & (df["L"] == 32)]
    df021 = df021.sort_values("shaping_rate")
    shaping_rates = 1 -df021["shaping_rate"]
    res = [nom_kl(r, p, p02) for r,p in zip(shaping_rates, df021["obtained_p"])]
    plt.plot(shaping_rates, res, label="p=0.24 L=32", marker="*", color="red")
    df021 = df[(df["target_p"] == p03) & (df["L"] == 1)]
    df021 = df021.sort_values("shaping_rate")
    shaping_rates = 1 -df021["shaping_rate"]
    res = [nom_kl(r, p, p03) for r,p in zip(shaping_rates, df021["obtained_p"])]
    plt.plot(shaping_rates, res, label="p=0.35 L=1", marker="*", color="cornflowerblue")
    df021 = df[(df["target_p"] == p03) & (df["L"] == 32)]
    df021 = df021.sort_values("shaping_rate")
    shaping_rates = 1 -df021["shaping_rate"]
    res = [nom_kl(r, p, p03) for r,p in zip(shaping_rates, df021["obtained_p"])]
    plt.plot(shaping_rates, res, label="p=0.24 L=32", marker="*", color="blue")

    h22 = h2(p02)
    h23 = h2(p03)
    plt.axvline(h22, linestyle="--", label="h2(0.24)", color="darkred")
    plt.axvline(h23, linestyle="--", label="h2(0.35)", color="darkblue")


    plt.grid(True, which="both")
    plt.xlabel(rf"$R_{{S}}$", fontsize=fs)
    #plt.ylabel(rf"$\frac{{1}}{{N}}\mathbb{{D}}(\tilde{{X}}_1^{{N}} || X_1^{{N}})$ [log]", fontsize=fs)
    plt.ylabel("Normalized Infromational Divergence [log]", fontsize=fs)
    plt.yscale("log")
    plt.legend(fontsize=fs)
    plt.show()

def plotcrcrm_comp():
    fs = 20
    #plot_for_id(417, color="#d00000ff", mark="X", leg=f"7-CRC-PC SCL-32&8")
    #plot_for_id(416, color="#0000d5ff", mark="X", leg=f"RowMerge-PC SCL-32&8")
    #plot_for_id(418, color="#d00000ff", mark="^", leg=f"7-CRC-PC SCL-32&32")
    #plot_for_id(138, color="#0000d5ff", mark="^", leg=f"RowMerge-PC SCL-32&32")
    plot_for_id(420, color="#0000d5ff", mark="^", leg=f"RowMerge-PC SCL-32&32")
    plot_for_id(421, color="#d00000ff", mark="^", leg=f"7-CRC-PC SCL-32&32")
    
    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2*fs)
    plt.xlabel("SNR [dB]", fontsize=2*fs)
    plt.tick_params(axis="both",labelsize=fs)
    #plt.xlim(2,7)
    #mplcursors.cursor(hover=True)
    plt.yscale("log")
    plt.legend(fontsize=fs)
    plt.show()

def plotcrc():
    fs = 20
    plot_for_id(263, color="#d0a000ff", mark="^", leg=r"3-CRC&3-CRC-PC SCL-1&32")
    plot_for_id(267, color="#007879ff", mark="^", leg=r"3-CRC&3-CRC-PC SCL-1&32 ")
    plot_for_id(271, color="#0000d5ff", mark="^", leg=r"3-CRC&3-CRC-PC SCL-8&8")
    plot_for_id(279, color="#d00000ff", mark="^", leg=r"5-CRC&5-CRC-PC SCL-32&32")

    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2*fs)
    plt.xlabel("SNR [dB]", fontsize=2*fs)
    plt.tick_params(axis="both",labelsize=fs)
    #plt.xlim(2,7)
    #mplcursors.cursor(hover=True)
    plt.yscale("log")
    plt.legend(fontsize=fs)
    plt.show()

def plotrm():
    fs = 20
    plot_for_id(316, color="#007879ff", mark="^", leg=r"1-RowMerge&8-RowMerge-CRC-PC SCL-1&32 ")
    plot_for_id(325, color="#d0a000ff", mark="^", leg=r"10-RowMerge&7-RowMerge-CRC-PC SCL-8&8")
    plot_for_id(312, color="#0000d5ff", mark="^", leg=r"1-RowMerge&5-CRC-PC SCL-8&32")
    plot_for_id(315, color="#d00000ff", mark="^", leg=r"1-RowMerge&5-CRC-PC SCL-32&32")

    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2*fs)
    plt.xlabel("SNR [dB]", fontsize=2*fs)
    plt.tick_params(axis="both",labelsize=fs)
    #plt.xlim(2,7)
    #mplcursors.cursor(hover=True)
    plt.yscale("log")
    plt.legend(fontsize=fs)
    plt.show()

def plotpc():
    fs = 20
    plot_for_id(370, color="#d00000ff", mark="^", leg=r"6-ParCheck&11-RowMerge-PC SCL-8&8")
    plot_for_id(352, color="#0000d5ff", mark="^", leg=r"4-ParCheck&5-CRC-PC 8&32")
    plot_for_id(360, color="#007879ff", mark="^", leg=r"4-ParCheck&5-CRC-PC SCL-32&32")
    plot_for_id(376, color="#d0a000ff", mark="^", leg=r"4-ParCheck&10-RowMerge-PC SCL-32&32")

    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2*fs)
    plt.xlabel("SNR [dB]", fontsize=2*fs)
    plt.tick_params(axis="both",labelsize=fs)
    #plt.xlim(2,7)
    #mplcursors.cursor(hover=True)
    plt.yscale("log")
    plt.legend(fontsize=fs)
    plt.show()

def plotfinal():
    fs = 20
    plot_for_id(279,color="#009564", mark="^", leg=r"5-CRC&5-CRC-PC SCL-32&32 ($R_M = 0.79$)")
    plot_for_id(315,color="#0000b6c4", mark="^", leg=r"1-RowMerge&5-CRC-PC SCL-32&32 ($R_M = 0.72$)")
    plot_for_id(225,color="#7e00ffc4", mark="^", leg=r"5-CRC-PC SCL-32&32 ($R_M = 0.79$)")
    plot_for_id(376, color="#fc8500", mark="^", leg=r"4-ParCheck&10-RowMerge-PC SCL-32&32 ($R_M = 0.7$)")
    plot_for_id(360, color="#ff0059", mark="^", leg=r"4-ParCheck&5-CRC-PC SCL-32&32 ($R_M = 0.7$)")

    plt.grid(True, which="both")
    plt.ylabel("FER", fontsize=2*fs)
    plt.xlabel("SNR [dB]", fontsize=2*fs)
    plt.tick_params(axis="both",labelsize=fs)
    #plt.xlim(2,7)
    mplcursors.cursor(hover=True)
    plt.yscale("log")
    plt.legend(fontsize=fs)
    plt.show()


#CRC + RM : 141
#CRC + CRC : 151
#CRC + PC : 152
#CRC + none : 117   (-119)
#CRC + no-shaping : 

#plot_rate_loss(range(88, 104))
#plot_smt()
#plot_rate_loss2()
#plot_for_id(106, toplot="fer")
#plot_for_ids(range(128, 138))
#139 baleze !! 
#plot_for_ids([146,147,148, 117])
#pour demain, plot à partir de 143
#plot_for_ids([117,141,151, 152])
#plot_for_ids([138, 139])
#plot_shaping_perf()
#plot_shaping_perf()
#plot_normalized_info_div()
#plot_for_shaping_pc()
#plot_for_ids([138, 410, 411, 412, 413, 414, 415, 416, 417, 418])
#plot_for_ids([420, 421])
#plotcrcrm_comp()
#plot_for_ids([225, 242])
plotfinal()
#plotrm()
#plotcrc()
#plotpc()
#plot_for_ids([i for i in range(190,422)])
#plot_shaping_perf_for_ids(range(190, 422))

#plot_for_ids([360, 353, 376, 361, 380, 352, 225, 357, 372, 365, 373, 315, 312, 309, 368, 374, 382])
#plot_for_ids([296, 279, 225, 360, 353, 365, 376, 315, 382])
#plot_shaping_perf_for_ids([279, 315, 225, 376, 360])
#plot_for_ids([315, 376, 360, 279, 225])
#plot_shaping_perf()
#plot_for_ids2([279, 315, 225, 376, 360])
#get_lower_bound_on_l1_dist_sm3()
