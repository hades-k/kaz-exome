import re
import numpy as np
import matplotlib.pyplot as plt

text = """Distance method: Pairwise difference
                     1         2         3         4         5         6         7         8         9        10        11        12        13
           1   0.00000
           2   0.00467   0.00000
           3   0.00396   0.00375   0.00000
           4   0.00243   0.00481   0.00411   0.00000
           5   0.04099   0.03602   0.02996   0.04323   0.00000
           6   0.04023   0.03676   0.03127   0.04331   0.00706   0.00000
           7   0.04081   0.03598   0.03003   0.04267   0.00265   0.01083   0.00000
           8   0.03994   0.03450   0.02874   0.04149   0.00420   0.01253   0.00164   0.00000
           9   0.05802   0.07578   0.07293   0.06994   0.10824   0.10105   0.10743   0.10747   0.00000
          10   0.05579   0.07296   0.07000   0.06759   0.10391   0.09627   0.10340   0.10347   0.00951   0.00000
          11   0.05813   0.07474   0.07218   0.06943   0.10645   0.09852   0.10560   0.10598   0.01800   0.00699   0.00000
          12   0.05301   0.07034   0.06751   0.06474   0.10292   0.09590   0.10234   0.10229   0.00192   0.00659   0.01488   0.00000
          13   0.02208   0.02930   0.02530   0.02899   0.03920   0.03523   0.03978   0.03975   0.03600   0.02648   0.02833   0.03142   0.00000
"""

def parse_lower_triangular(text):
    data_lines = []
    for line in text.splitlines():
        if re.match(r'^\s*\d+\s', line):
            data_lines.append(line.strip())

    row_indices = []
    row_values = []
    for line in data_lines:
        parts = re.split(r'\s+', line)
        idx = int(parts[0])
        vals = list(map(float, parts[1:]))
        row_indices.append(idx)
        row_values.append(vals)

    n = max(row_indices)
    M = np.zeros((n, n), dtype=float)

    for idx, vals in zip(row_indices, row_values):
        i = idx - 1
        for j, v in enumerate(vals):
            M[i, j] = v
            M[j, i] = v

    np.fill_diagonal(M, 0.0)
    labels = [str(i) for i in range(1, n + 1)]
    return M, labels

M, labels = parse_lower_triangular(text)

custom_labels = ['Bengali', 'Gujarati', 'Paunjabi',
                 'Tamil', 'British', 'Finnish',
                 'Iberian', 'Toscani', 'Dai Chinese',
                 'Han Chinese', 'Japanese',
                 'Kinh Vietnamese', 'Kazakh']

fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(M, cmap='YlGn', interpolation='nearest')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
plt.rcParams.update({'mathtext.default':  'regular' })
cbar.set_label('$F_{ST}$')

ax.set_xticks(np.arange(len(labels)))
ax.set_yticks(np.arange(len(labels)))
ax.set_xticklabels(custom_labels, rotation=90)
ax.set_yticklabels(custom_labels)
plt.rcParams.update({'mathtext.default':  'regular' })
ax.set_title('Matrix of pairwise $F_{ST}$')
ax.set_aspect('equal')
fig.tight_layout()
plt.show()
