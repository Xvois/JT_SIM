import pandas as pd
import matplotlib.pyplot as plt

# Matplotlib settings
plt.rcParams['figure.figsize'] = [6, 4]  # Adjust the width to 8 inches
plt.rcParams['figure.dpi'] = 100

# Interior ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# Make the font size larger
plt.rcParams.update({'font.size': 12})

df = pd.read_csv('ChamberTemps[IDEAL].csv')

# A function that when a NA is spotted duplicates the previous value
def fill_na(series):
    series = series.copy()  # Create a copy to avoid modifying the original DataFrame
    for i in range(1, len(series)):
        if pd.isna(series[i]):
            series.loc[i] = series[i - 1]
    return series

delta_T = df['RT'] - df['LT']
delta_P = df['RP'] - df['LP']

length = len(delta_T)



LT_moving_avg = fill_na(df['LT']).rolling(window=(length // 10)).mean()
RT_moving_avg = fill_na(df['RT']).rolling(window=(length // 10)).mean()
ET_moving_avg = fill_na(df['ET']).rolling(window=(length // 10)).mean()

plt.plot(df['time'], df['LT'], color='grey', alpha=0.1)
plt.plot(df['time'], LT_moving_avg, label='Left Temp Moving Avg', color='blue')
plt.plot(df['time'], df['RT'], color='grey', alpha=0.1)
plt.plot(df['time'], RT_moving_avg, label='Right Temp Moving Avg', color='red')
plt.plot(df['time'], df['ET'], color='grey', alpha=0.1)
plt.plot(df['time'], ET_moving_avg, label='Ensemble Temp Moving Avg', color='green')
plt.xlabel('time')
plt.ylabel('Temperature')
plt.ylim(min(df['ET']), max(df['ET']))
plt.legend()
plt.savefig('T_sys.png')
plt.show()

delta_T_mov_avg = RT_moving_avg - LT_moving_avg

plt.plot(df['time'], delta_T, color='grey', alpha=0.1)
plt.plot(df['time'], delta_T_mov_avg, label='Temperature Difference', color='purple')
plt.xlabel('time')
plt.ylim(-5,5)
plt.ylabel('Temperature Difference')
plt.savefig('T_diff.png')
plt.show()

RP_moving_avg = fill_na(df['RP']).rolling(window=(length // 10)).mean()
LP_moving_avg = fill_na(df['LP']).rolling(window=(length // 10)).mean()

delta_P_mov_avg = RP_moving_avg - LP_moving_avg

plt.plot(df['time'], delta_P_mov_avg, label='Left Pressure Moving Avg', color='blue')
plt.plot(df['time'], delta_P, color='grey', alpha=0.1)
plt.xlabel('time')
plt.ylabel('Pressure Difference')
plt.savefig('P_diff.png')
plt.show()

fig, ax1 = plt.subplots()

color = 'tab:blue'
ax1.set_ylabel('Pressure Difference  / nNm$^{-1}$', color=color)
ax1.plot(df['time'] * 1e+6, delta_P_mov_avg * 1e+9, label='Pressure Difference', color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xlabel('Time / $\mu$s')

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel('Temperature Difference  / K', color=color)
ax2.plot(df['time'] * 1e+6, delta_T_mov_avg, label='Temperature Difference', color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(-5, 25)
ax2.set_xlabel('Time / $\mu$s')

fig.tight_layout()
plt.savefig('P_T_diff.png')
plt.show()

# Calculate the final JT value
final_JT = (df['RT'].iloc[-1] - df['LT'].iloc[-1]) / (df['RP'].iloc[-1] - df['LP'].iloc[-1])
print(f"Final JT value: {final_JT}")
