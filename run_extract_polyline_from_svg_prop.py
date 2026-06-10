import numpy as np
import matplotlib.pyplot as plt
from svgpathtools import svg2paths
from my_modules import cubic_bezier

fn = './data_files/propeller_diagram_clean1.svg'

paths, attributes = svg2paths(fn)
n_paths = len(paths)
n_attributes = len(attributes)
print('n paths = ',n_paths)

profiles = []
for j in range(n_paths):
    attr = attributes[j]
    color = attr.get('stroke','1')
    if color=='rgb(100%, 0%, 100%)':
        iseg = 0
        profile = []
        for segment in paths[j]:
            if type(segment).__name__=='Line':
                line = [[np.imag(segment.start),np.imag(segment.end)],\
                        [np.real(segment.start),np.real(segment.end)]]
                profile.append(line)
                
            elif type(segment).__name__=='CubicBezier':
                cp = [
                    (np.imag(segment.start),np.real(segment.start)),
                    (np.imag(segment.control1),np.real(segment.control1)),
                    (np.imag(segment.control2),np.real(segment.control2)),
                    (np.imag(segment.end),np.real(segment.end))
                ]
                profile.append(cp)
        if len(profile)>0:
            profiles.append(profile)

t = np.linspace(0, 1, 91)
sections = []
nprof = len(profiles)
print('n profiles = ',nprof)
for i in range(nprof):
    section = np.zeros([0,2],float)
    profile = profiles[i]
    nseg = len(profile)
    for j in range(nseg):
        seg = profile[j]
        nline = len(seg)

        if nline==2:
            add_seg = np.array(seg).T
        elif nline==4:
            add_seg = np.array([cubic_bezier(ti, np.array(seg[0]), np.array(seg[1]), np.array(seg[2]), np.array(seg[3])) for ti in t])
        else:
            add_seg = np.zeros([0,2],float)

        if section.shape[0]==0:
            section = np.append(section,add_seg,axis=0)
        elif np.array_equal(section[-1,:],add_seg[0,:]):
            section = np.append(section,add_seg[1:,:],axis=0)
        else:
            section = np.append(section,add_seg,axis=0)

    sections.append(section)

sections = sections[::-1]

for i in range(12):
    fn = './data_files/sections_propeller_sec{:02d}.dat'.format(i)
    fmt='%.18e'
    np.savetxt(fn,sections[i],fmt=fmt,delimiter='\t')

if 0:
    rad = np.array([72.5,87.0,101.5,116.0,145.0,174.0,203.0,232.0,261.0,275.5,282.8,287.1])
    fig, ax = plt.subplots()
    for i in range(12):
        sec = sections[i]
        ax.plot(sec[:,0],sec[:,1]+rad[i])

    plt.axis('equal')
    plt.show()