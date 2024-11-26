import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('ChamberTemps.csv')

delta = df['RT'] - df['LT']

length = len(delta)

LT_moving_avg = df['LT'].rolling(window=(length // 20)).mean()
RT_moving_avg = df['RT'].rolling(window=(length // 20)).mean()

plt.plot(df['Iteration'], df['LT'], color='grey', alpha=0.1)
plt.plot(df['Iteration'], LT_moving_avg, label='Left Temp Moving Avg', color='blue')
plt.plot(df['Iteration'], df['RT'], color='grey', alpha=0.1)
plt.plot(df['Iteration'], RT_moving_avg, label='Right Temp Moving Avg', color='red')
plt.plot(df['Iteration'], df['ET'], label='Ensemble Temp [SHOULD BE CONST]', color='green')
plt.xlabel('Iteration')
plt.ylabel('Temperature')
plt.legend()
plt.savefig('T_sys.png')
plt.show()

delta_mov_avg = RT_moving_avg - LT_moving_avg

plt.plot(df['Iteration'], delta, color='grey', alpha=0.1)
plt.plot(df['Iteration'], delta_mov_avg, label='Temperature Difference', color='purple')
plt.xlabel('Iteration')
plt.ylabel('Temperature Difference')
plt.savefig('T_diff.png')
plt.show()
