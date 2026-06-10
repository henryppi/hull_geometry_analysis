import numpy as np
import matplotlib as plt
# points = np.array(points)
# points = np.unique(points,axis=0)

# flip direction
# points[:,0] *= -1

# for cubic in cubics:
#     ax.add_patch(cubic)
# ax.plot(points[:,0],points[:,1],'.-k')

# rad = np.array([72.5,87.0,101.5,116.0,145.0,174.0,203.0,232.0,261.0,275.5,282.8,287.1])
# fig, ax = plt.subplots()
# for i in range(12):
#     sec = sections[i]
#     ax.plot(sec[:,0],sec[:,1]+rad[i])

# plt.axis('equal')
# plt.show()

rad = np.array([72.5,87.0,101.5,116.0,145.0,174.0,203.0,232.0,261.0,275.5,282.8,287.1])
tan_lead = np.array([36,54,71,86,111,124,124,105,59,18,-13,-41])
tan_tail = np.array([35,46,57,67,87,107,128,149,165,165,145,103])