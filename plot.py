import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('ChamberTemps.csv')

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
print(LT_moving_avg)
RT_moving_avg = fill_na(df['RT']).rolling(window=(length // 10)).mean()

plt.plot(df['time'], df['LT'], color='grey', alpha=0.1)
plt.plot(df['time'], LT_moving_avg, label='Left Temp Moving Avg', color='blue')
plt.plot(df['time'], df['RT'], color='grey', alpha=0.1)
plt.plot(df['time'], RT_moving_avg, label='Right Temp Moving Avg', color='red')
plt.plot(df['time'], df['ET'], label='Ensemble Temp', color='green')
plt.xlabel('time')
plt.ylabel('Temperature')
plt.ylim(0, 400)
plt.legend()
plt.savefig('T_sys.png')
plt.show()

delta_T_mov_avg = RT_moving_avg - LT_moving_avg

plt.plot(df['time'], delta_T, color='grey', alpha=0.1)
plt.plot(df['time'], delta_T_mov_avg, label='Temperature Difference', color='purple')
plt.xlabel('time')
plt.ylim(-50,50)
plt.ylabel('Temperature Difference')
plt.savefig('T_diff.png')
plt.show()

RP_moving_avg = fill_na(df['RP']).rolling(window=(length // 20)).mean()
LP_moving_avg = fill_na(df['LP']).rolling(window=(length // 20)).mean()

delta_P_mov_avg = RP_moving_avg - LP_moving_avg

plt.plot(df['time'], delta_P_mov_avg, label='Left Pressure Moving Avg', color='blue')
plt.plot(df['time'], delta_P, color='grey', alpha=0.1)
plt.xlabel('time')
plt.ylabel('Pressure Difference')
plt.savefig('P_diff.png')
plt.show()

# Remove large outliers
JT = delta_T_mov_avg / delta_P_mov_avg

# Get last JT value and print
JT_last = JT.iloc[-1]
print(f'Last JT value: {JT_last}')

plt.show()
