import numpy as np
import matplotlib.pyplot as plt
import csv

with open('./data.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

# print(data[2][1])
i = 0

time = []
nis_laser = []
thresh_laser = []
nis_radar = []
thresh_radar = []

headings = data[0]
data = data[1:]
# print(data[0])
# print(data[0][0])

for category in data:
    time.append((int(category[0]) - int(data[0][0])) / 1000000)
    nis_laser.append(float(category[1]))
    thresh_laser.append(float(category[2]))
    nis_radar.append(float(category[3]))
    thresh_radar.append(float(category[4]))

# print(headings[1])

plt.subplot(2,1,1)
plt.plot(time, nis_laser, label = headings[1])
plt.plot(time, thresh_laser, label = headings[2])
plt.xlabel(headings[0])
plt.ylabel('NIS Laser')
plt.yticks(np.arange(0, 21, step = 5.0))
plt.legend()

plt.subplot(2,1,2)
plt.plot(time, nis_radar, label = headings[3])
plt.plot(time, thresh_radar, label = headings[4])
plt.xlabel(headings[0])
plt.ylabel('NIS Radar')
plt.yticks(np.arange(0, 21, step = 5.0))
plt.legend()

plt.show()
