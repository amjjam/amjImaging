import sys
sys.path.append("../src/")
import image

img=image.read("../dat/test1.dat")

print(img)

import matplotlib.pyplot as plt
plt.imshow(img)

plt.show()
