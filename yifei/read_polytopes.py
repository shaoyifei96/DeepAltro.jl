import numpy as np
path = '/home/yifei/Documents/optimal_ctrl/Altro.jl/yifei/data/h_poly1.txt'
with open(path, 'r') as file:
    data = file.read().split('\n\n')
    # print(data)
for data_i in data:
    if len(data_i) > 1:
        polytope = np.matrix(data_i.replace(',',' ').replace('\n',';'))
        print(polytope.shape)

