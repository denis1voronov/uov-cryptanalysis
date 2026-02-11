import math
import numpy as np
import matplotlib as mpl
mpl.rcParams.update({"mathtext.fontset": "cm", "axes.unicode_minus": False,})
import matplotlib.pyplot as plt
from cryptographic_estimators.MQEstimator import MQEstimator

uov1 = dict(name=r"$\mathtt{uov\!-\!Ip}$",  v=68, o=44, m=44, q=256, b=8,
            color="green", marker="o")
uov2 = dict(name=r"$\mathtt{uov\!-\!Is}$",  v=96, o=64, m=64, q=16, b=4,
            color="orange", marker="o")
qruov1 = dict(name=r"$\mathtt{qruov\!-\!Ia}$",  v=55, o=20, m=60, q=31**3, b=15,
              color="purple", marker="o")
qruov2 = dict(name=r"$\mathtt{qruov\!-\!Ib}$",  v=60, o=7, m=70, q=31**10, b=50,
              color="red", marker="o")
mayo1 = dict(name=r" $\mathsf{MAYO}_1$",  v=78, o=8, m=78, q=16, b=4,
             color="teal", marker="o")
mayo2 = dict(name=r" $\mathsf{MAYO}_2$",  v=64, o=17, m=64, q=16, b=4,
             color="olive", marker="o")
number = 100

def mq_bit_complexity(var: int, m: int, q: int):
    instance = MQEstimator(n=var, m=m, q=q).booleansolve_fxl
    return instance.time_complexity(k=instance.k(), variant="las_vegas")

def cache(m: int, q: int):
    complexities = {}
    for k in range(1, m + 1): # we only need complexities of overdetermined mq
        complexities[k] = mq_bit_complexity(k, m, q)
    return complexities

def estimate(w1: int, *, m: int, cache: dict[int, float],
             ers: np.ndarray, rtrn: bool = False):
    bits = 0.0
    bad = 0
    values = [] if rtrn else None

    for i in range(number):
        erasures = ers[i]
        column_sums = erasures.sum(axis=0)         # count erasures per column
        best_column = int(np.argmin(column_sums))  # choose column with min erasures
        omega = erasures[:, best_column]           # erasures vector

        erased = (omega > 0) # choose erased variables
        w = int(erased.sum()) # compute w = |W|
        if w == 0:
            if w1 == 0:
                attack_bits = 0.0 # no enumeration
                bits += attack_bits
                if rtrn:
                    values.append(attack_bits)
            else:
                bad += 1
            continue

        if w1 > w:
            bad += 1
            continue

        var = w - w1 # choose mq variables
        if var < 0 or var > m:
            bad += 1
            continue

        if var == 0:
            mq_bits = 0.0
        else:
            mq_bits = cache[var]

        w_erased = omega[erased]
        w_sorted = np.sort(w_erased)

        if w1 == 0:
            enum_bits = 0
        else:
            enum_bits = int(w_sorted[:w1].sum())

        attack_bits = enum_bits + mq_bits
        bits += attack_bits
        if rtrn:
            values.append(attack_bits)

    good = number - bad
    if good == 0:
        return math.inf, (np.array([], dtype=float) if rtrn else None)

    avg_bits = bits / good
    if rtrn:
        return avg_bits, np.array(values, dtype=float)
    return avg_bits, None

def curve(params: dict, p_range: np.ndarray, rng: np.random.Generator):
    v, m, o, q, b = (params["v"], params["m"],
    params["o"], params["q"], params["b"])
    mq = cache(m, q)
    average_complexities: list[float] = []
    best_parameters: list[int | None] = []
    all_complexities: list[np.ndarray] = []

    print(f"[{params['name']}] p (average bit complexity, w~)")
    for p in p_range: # generates bit erasures in v*o matrix "number" times
        erasures = rng.binomial(n=b, p=p, size=(number, v, o))
        best_complexity = math.inf
        best_parameter = None

        for w1 in range(0, v + 1):
            bits, _ = estimate(w1, m=m, cache=mq, ers=erasures, rtrn=False)
            if bits < best_complexity:
                best_complexity = bits
                best_parameter = w1

        if best_parameter is None or not math.isfinite(best_complexity):
            values = np.array([], dtype=float)
        else:
            _, values = estimate(best_parameter, m=m, cache=mq,
            ers=erasures, rtrn=True)

        average_complexities.append(best_complexity)
        best_parameters.append(best_parameter)
        all_complexities.append(values)

        print(f"{p:.2f} ({best_complexity:.2f}, {best_parameter})")

    return average_complexities, best_parameters, all_complexities

def main():
    rng = np.random.default_rng(0)
    p_range = np.arange(0.0, 1.001, 0.01)
    # uov1_curve, uov1_w, uov1_values = curve(uov1, p_range, rng)
    # uov2_curve, uov2_w, uov2_values = curve(uov2, p_range, rng)
    # qruov1_curve, qruov1_w, qruov1_values = curve(qruov1, p_range, rng)
    # qruov2_curve, qruov2_w, qruov2_values = curve(qruov2, p_range, rng)
    mayo1_curve, mayo1_w, mayo1_values = curve(mayo1, p_range, rng)
    mayo2_curve, mayo2_w, mayo2_values = curve(mayo2, p_range, rng)
    plt.figure(figsize=(4.4, 3.3))

    def plot_scheme_dots(p_range, scheme_values, scheme): # plots all dots
        for p, vals in zip(p_range, scheme_values):
            if vals.size:
                plt.scatter(np.full(vals.shape, p, dtype=float), vals, s=4,
                marker=scheme["marker"], color=scheme["color"],
                alpha=0.05, linewidths=0, zorder=1)

    # plot_scheme_dots(p_range, uov1_values, uov1)
    # plot_scheme_dots(p_range, uov2_values, uov2)
    # plot_scheme_dots(p_range, qruov1_values, qruov1)
    # plot_scheme_dots(p_range, qruov2_values, qruov2)
    plot_scheme_dots(p_range, mayo1_values, mayo1)
    plot_scheme_dots(p_range, mayo2_values, mayo2)

    def plot_scheme_avg(p_range, scheme_curve, scheme): # plots average dots
        plt.scatter(p_range, scheme_curve, marker=scheme["marker"], s=4,
                    color=scheme["color"], label=scheme["name"])

    # plot_scheme_avg(p_range, uov1_curve, uov1)
    # plot_scheme_avg(p_range, uov2_curve, uov2)
    # plot_scheme_avg(p_range, qruov1_curve, qruov1)
    # plot_scheme_avg(p_range, qruov2_curve, qruov2)
    plot_scheme_avg(p_range, mayo1_curve, mayo1)
    plot_scheme_avg(p_range, mayo2_curve, mayo2)

    plt.xlim(0, max(p_range))
    plt.ylim(0, 350)
    plt.margins(x=0, y=0)
    plt.axhline(143, color="black", linestyle="--", linewidth=1)
    plt.xlabel("bit erasure probability $p$", fontsize=12)
    plt.ylabel("bit complexity", fontsize=12)
    plt.grid(True, which="major", linewidth=1.2, alpha=0.7)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax = plt.gca()
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: rf"${x:g}$"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: rf"${x:g}$"))
    plt.tight_layout()
    plt.savefig("attack.pdf", bbox_inches="tight")
    plt.show()

    # uov1_w_arr = np.array(uov1_w, dtype=float)
    # uov2_w_arr = np.array(uov2_w, dtype=float)
    # qruov1_w_arr = np.array(qruov1_w, dtype=float)
    # qruov2_w_arr = np.array(qruov2_w, dtype=float)
    mayo1_w_arr = np.array(mayo1_w, dtype=float)
    mayo2_w_arr = np.array(mayo2_w, dtype=float)
    plt.figure(figsize=(4.4, 3.3))

    def plot_scheme_w(p_range, scheme_arr, scheme): # plots best w~
        plt.scatter(p_range, scheme_arr,
        marker="o", s=4, color=scheme["color"], label=scheme["name"])

    # plot_scheme_w(p_range, uov1_w_arr, uov1)
    # plot_scheme_w(p_range, uov2_w_arr, uov2)
    # plot_scheme_w(p_range, qruov1_w_arr, qruov1)
    # plot_scheme_w(p_range, qruov2_w_arr, qruov2)
    plot_scheme_w(p_range, mayo1_w_arr, mayo1)
    plot_scheme_w(p_range, mayo2_w_arr, mayo2)

    plt.xlim(0, max(p_range))
    plt.ylim(0, 60)
    plt.margins(x=0, y=0)
    plt.xlabel("bit erasure probability $p$", fontsize=12)
    plt.ylabel(r"best $\tilde{w}$ parameter", fontsize=12)
    plt.grid(True, which="major", linewidth=1.2, alpha=0.7)
    plt.legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax = plt.gca()
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: rf"${x:g}$"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: rf"${x:g}$"))
    plt.tight_layout()
    plt.savefig("attack_bestparams.pdf", bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    main()