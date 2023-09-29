import numpy as np


def adaptive_walk(T, local_peaks=None, max_steps=5, conv_crit=10**-5, sample_size=10, starting_set=None):
    n = T.shape[0]
    paths = {}
    counts = []
    
    for s in starting_set:
        i = 1
        while i <= sample_size:
            path = [s]
            node = s
            j = 1
            while j < max_steps and not local_peaks[node]:
                node = np.random.choice(T[node].indices, p=T[node].data, size=1)[0]
                path.append(node)
                j += 1
            path = tuple(path)
            if path in paths:
                counts[paths[path]] += 1
            else:
                paths[path] = len(paths)
                counts.append(1)
            i += 1
        
    return paths, counts