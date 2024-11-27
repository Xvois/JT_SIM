import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('ChamberTemps.csv')

delta_T = df['RT'] - df['LT']
delta_P = df['RP'] - df['LP']

length = len(delta_T)

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

delta_T_mov_avg = RT_moving_avg - LT_moving_avg

plt.plot(df['Iteration'], delta_T, color='grey', alpha=0.1)
plt.plot(df['Iteration'], delta_T_mov_avg, label='Temperature Difference', color='purple')
plt.xlabel('Iteration')
plt.ylabel('Temperature Difference')
plt.ylim(-20, 20)
plt.savefig('T_diff.png')
plt.show()

RP_moving_avg = df['RP'].rolling(window=(length // 20)).mean()
LP_moving_avg = df['LP'].rolling(window=(length // 20)).mean()

delta_P_mov_avg = RP_moving_avg - LP_moving_avg

plt.plot(df['Iteration'], delta_P_mov_avg, label='Left Pressure Moving Avg', color='blue')
plt.plot(df['Iteration'], delta_P, color='grey', alpha=0.1)
plt.xlabel('Iteration')
plt.ylabel('Pressure Difference')
plt.savefig('P_diff.png')
plt.show()

plt.plot(df['Iteration'], delta_T_mov_avg / delta_P_mov_avg, label='JT', color='purple')
plt.xlabel('Iteration')
plt.ylabel('Joule-Thomson Coefficient')

plt.savefig('JT.png')
plt.show()
