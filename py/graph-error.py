import matplotlib.pyplot as plt
import pandas as pd

far = pd.read_csv("./WORKSPACE/calc/3-11-2025_8h7m6s_gmm128_TV100/scoreFAR")
frr = pd.read_csv("./WORKSPACE/calc/3-11-2025_8h7m6s_gmm128_TV100/scoreFRR")

plt.hist(far, 40)
plt.hist(frr, 40)
plt.show()
