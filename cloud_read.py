from scipy.io import FortranFile
import matplotlib.pyplot as plt
from constants import *
import numpy as np

f = FortranFile("outputdata1/nube3101.sal", "r")
data = f.read_reals(np.float32)
U1 = data[: biased_nx1 * biased_nx1 * biased_nz1].reshape(
    (biased_nx1, biased_nx1, biased_nz1), order="F"
)
U2 = data[
    biased_nx1 * biased_nx1 * biased_nz1: 2 * biased_nx1 * biased_nx1 *
    biased_nz1].reshape((biased_nx1, biased_nx1, biased_nz1), order="C")
plt.imshow(U1[:, :, 10])
plt.savefig("img.png")
print(len(data) / (biased_nx1 * biased_nx1 * biased_nz1))
# .reshape((57, 57, 49), order="F")
