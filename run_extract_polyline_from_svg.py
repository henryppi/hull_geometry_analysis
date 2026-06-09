# extract images
# pdfimages -all 16229.pdf rudder_ 

# extract svg  
# pdftocairo -svg 16229.pdf output.svg

import numpy as np
import matplotlib.pyplot as plt 

from svgpathtools import svg2paths
fn = './data_files/output_clean_profil2.svg'
# paths is a list of Path objects; attributes is a list of dictionaries for each tag
paths, attributes = svg2paths(fn)
n_paths = len(paths)
n_attributes = len(attributes)
print('n paths = ',n_paths)
print('n attributes = ',n_attributes)
if 0:
    # To get the raw 'd' attribute string for each path
    for attr in attributes:
        print(attr['d'])

center = []
points = []
for j in range(n_paths):
    # To access the geometry (Line, CubicBezier, etc.) of the first path
    a = []
    b = []
 
    iseg = 0
    for segment in paths[j]:
        attr = attributes[j]
        lw = float(attr.get('stroke-width', '1'))
        
        if 0:
            print('#path = ',j,'\t #seg = ',iseg, '\t lw = ',lw)
            print(f"Start: {segment.start}, End: {segment.end}")

        if lw>1.0:
            points.append([np.real(segment.start),np.imag(segment.start)])
            points.append([np.real(segment.end),np.imag(segment.end)])
        if lw<1.0:
            center.append([np.real(segment.start),np.imag(segment.start)])
            center.append([np.real(segment.end),np.imag(segment.end)])
        # print(center)
        iseg+=1


center = np.array(center)
xmax = center[np.argmax(center[:,0]),:]
xmin = center[np.argmin(center[:,0]),:]
ymax = center[np.argmax(center[:,1]),:]
ymin = center[np.argmin(center[:,1]),:]
cc = np.array([0.5*(ymax[0]+ymin[0]),0.5*(xmax[1]+xmin[1])])

points = np.array(points)
points = np.unique(points,axis=0)

#align with csys
points[:,0] += -cc[0]
points[:,1] += -cc[1]

# flip direction
points[:,0] *= -1
plt.plot(points[:,0],points[:,1],'.k')

fn = './data_files/section_points_rudder.dat'
fmt='%.18e'
np.savetxt(fn,points,fmt=fmt,delimiter='\t')

# plt.axis('equal')
# plt.show()