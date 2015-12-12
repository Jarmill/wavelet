#piltest.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

img = mpimg.imread("..\\img\\dock.jpg")
fig1 = plt.imshow(img)

fig2 = plt.figure()
fig2 = plt.hist(img.ravel(), bins = 256, range=(0, 256))

from lift2d_test import *
iW = legall53_2d(np.copy(img[:, :, 1]), MAX_LAYER = 2)
W = np.dstack((iW, iW, iW))
fig3 = plt.imshow(W)
plt.show()

