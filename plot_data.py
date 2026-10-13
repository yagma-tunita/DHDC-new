import csv, os, sys
sys.stdout.reconfigure(encoding='utf-8')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'SimSun']
plt.rcParams['axes.unicode_minus'] = False

BASE = 'data'
OUT = 'png'
os.makedirs(OUT, exist_ok=True)

def read_csv(f):
    p = os.path.join(BASE, f)
    with open(p, 'rb') as fh:
        head = fh.read(2)
    enc = 'utf-16' if head in (b'\xff\xfe', b'\xfe\xff') else 'utf-8-sig'
    return list(csv.reader(open(p, encoding=enc)))

# ---- data_nrange.csv：四模型按服务器规模 N 对齐 ----
rows = read_csv('data_nrange.csv')
hdr = rows[0]
models = defaultdict(list)
for r in rows[1:]:
    if len(r) < 9 or not r[2].strip().isdigit():
        continue
    model = r[0]
    N = int(r[2]); degree = int(r[3]); diam = int(r[4])
    avghops = float(r[5]); loadfactor = float(r[6])
    maxthp = float(r[7]); avgthp = float(r[8])
    models[model].append(dict(N=N, degree=degree, diam=diam, avghops=avghops,
                              loadfactor=loadfactor, maxthp=maxthp, avgthp=avgthp))
for m in models:
    models[m].sort(key=lambda x: x['N'])

colors = {'DHDC': '#d62728', 'HSDC': '#1f77b4', 'CSDC': '#2ca02c', 'BCDC': '#ff7f0e'}
markers = {'DHDC': 'o', 'HSDC': 's', 'CSDC': '^', 'BCDC': 'v'}
order = ['DHDC', 'HSDC', 'CSDC', 'BCDC']

def plot_metric(key, ylabel, title, fname):
    fig, ax = plt.subplots(figsize=(7, 5))
    for m in order:
        xs = [p['N'] for p in models[m]]
        ys = [p[key] for p in models[m]]
        ax.plot(xs, ys, marker=markers[m], color=colors[m], label=m,
                linewidth=2, markersize=6)
    ax.set_xscale('log')
    ax.set_xlabel('服务器规模 N（对数坐标）', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=150)
    plt.close(fig)

plot_metric('degree', '端口数（Degree）', '单服务器端口数对比（按规模 N 对齐）', 'degree_vs_N.png')
plot_metric('loadfactor', '负载因子（loadFactor）', '全网拥塞对比（按规模 N 对齐）', 'loadfactor_vs_N.png')
plot_metric('avghops', '平均跳数（avgHops）', '平均延迟对比（按规模 N 对齐）', 'avghops_vs_N.png')
plot_metric('diam', '直径（Diameter）', '网络直径对比（按规模 N 对齐）', 'diameter_vs_N.png')
plot_metric('avgthp', '平均吞吐（Gbps）', '平均吞吐量对比（按规模 N 对齐）', 'avgthp_vs_N.png')

# ---- 吞吐-时间曲线（四模型对比）----
def plot_thp_compare(kind, fname, title):
    fig, ax = plt.subplots(figsize=(8, 5))
    for m in order:
        rows = read_csv('thp_%s_%s.csv' % (kind, m))
        xs = [float(r[0]) for r in rows[1:] if len(r) >= 2]
        ys = [float(r[1]) for r in rows[1:] if len(r) >= 2]
        ax.plot(xs, ys, color=colors[m], label=m, linewidth=1.5)
    ax.set_xlabel('时间（s）', fontsize=12)
    ax.set_ylabel('吞吐量（Gbps）', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=150)
    plt.close(fig)

plot_thp_compare('random', 'thp_random_curve.png', '四模型随机流量下的吞吐-时间曲线对比')
plot_thp_compare('alltoall', 'thp_alltoall_curve.png', '四模型全对全流量下的吞吐-时间曲线对比')

print('已生成图像：')
for f in sorted(os.listdir(OUT)):
    print(' ', os.path.join(OUT, f))
