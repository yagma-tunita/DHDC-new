import csv, os, sys
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

rows = read_csv('data_nrange.csv')
data = defaultdict(list)
for r in rows[1:]:
    if len(r) < 8 or not r[2].strip().isdigit():
        continue
    m = r[0]
    data[m].append(dict(N=int(r[2]), degree=int(r[3]), loadfactor=float(r[6])))
for m in data:
    data[m].sort(key=lambda x: x['N'])

colors = {'DHDC': '#d62728', 'HSDC': '#1f77b4', 'CSDC': '#2ca02c', 'BCDC': '#ff7f0e'}
markers = {'DHDC': 'o', 'HSDC': 's', 'CSDC': '^', 'BCDC': 'v'}
order = ['DHDC', 'HSDC', 'CSDC', 'BCDC']

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

for m in order:
    xs = [p['degree'] for p in data[m]]
    ys = [np.log10(p['N']) for p in data[m]]
    zs = [p['loadfactor'] for p in data[m]]
    lw = 3.2 if m == 'DHDC' else 2.0
    ax.plot(xs, ys, zs, color=colors[m], marker=markers[m], markersize=6,
            linewidth=lw, label=m, alpha=1.0 if m == 'DHDC' else 0.8)

ax.set_xlabel('端口数', fontsize=13, labelpad=10)
ax.set_ylabel('服务器规模 log10(N)', fontsize=13, labelpad=10)
ax.set_zlabel('拥塞 loadFactor', fontsize=13, labelpad=10)
ax.set_title('四种网络：端口数、服务器规模与拥塞的三维轨迹', fontsize=14)
ax.legend(fontsize=13, loc='upper left')
ax.view_init(elev=20, azim=-55)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'thp_3d_surface.png'), dpi=150)
plt.close(fig)
print('已生成 png/thp_3d_surface.png')
