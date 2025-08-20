import matplotlib.pyplot as plt
import csv

x = []
s11 = []
s21 = []
first = True

with open('postpro/coplanar_waveguide_driven_simulazione_lunga/port-S.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter=',')
    for row in lines:
        if first == False :
            x.append(float(row[0]))
            s11.append(float(row[1]))
            s21.append(float(row[3]))
        first = False

plt.plot(x, s11, color = 'g', linestyle = 'solid', marker = 'o',label = "S11")

plt.xticks(rotation = 25)
plt.xlabel('Frequency [GHz]')
plt.ylabel('S11 [dB]')
plt.title('S11 vs frequency', fontsize = 20)
plt.grid()
plt.legend()
plt.show()

plt.plot(x, s21, color = 'r', linestyle = 'solid', marker = 'o',label = "S21")

plt.xticks(rotation = 25)
plt.xlabel('Frequency [GHz]')
plt.ylabel('S21 [dB]')
plt.title('S21 vs frequency', fontsize = 20)
plt.grid()
plt.legend()
plt.show()
