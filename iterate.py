import numpy as np
import matplotlib.pyplot as plt





c, d = [], []

for i in xrange(1000):
    m, b = 0.6, -1
    x = np.arange(6, 14, 0.1)
    y = m*x+b
    noise = np.random.normal(0, 1, len(x))
    y = y + noise
    new_x = x[np.where(x > 7)]
    new_y = y[np.where(x > 7)]
    
    for num in xrange(10):
        print len(new_x), len(new_y)
        fit = np.polyfit(new_x, new_y, 1)
        m = fit[0]
        b = fit[1]
        scatter = np.std(new_y - m*new_x -b)
        cut = (3-b+scatter)/m
        new_x = x[np.where(x > cut)]
        new_y = y[np.where(x > cut)]
       
    c.append(m)
    d.append(b)

#print c
#print d
    
print np.mean(c), np.mean(d)


#print fit

#plt.plot(x, y, 'ko')
#plt.plot(x, m*x+b, 'r--')
#plt.hlines(3, 6, 14)
#plt.show()
