import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import viridis

mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.size'] = 1.5
mpl.rcParams['xtick.minor.width'] = 0.5
mpl.rcParams['ytick.major.size'] = 3
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.size'] = 1.5
mpl.rcParams['ytick.minor.width'] = 0.5

plt.rcParams['svg.fonttype'] = 'none'

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.serif'] = 'Arial'

cmap = plt.cm.viridis