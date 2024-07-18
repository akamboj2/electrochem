import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
# data = np.random.rand(60,60,60)
# print(data.shape)

# sheet = '4.4-2.5'
sheet = '4.56-2.5V'
data = pd.read_excel('IK-2-151C_SlowTimeSeries_000_300-1000.xlsx', sheet_name=sheet)


# Set the first column as the x axis
x_axis = data.iloc[:, 0]
print('x axis', x_axis.values)

# set the y axis as all the headers
y_axis = data.columns[1:]
print('y axis', y_axis)

# Remove the first column from the data frame
data = data.iloc[:, 1:]

# Transpose the data frame
data = data.transpose()

# # Create a heatmap with wave number on x axis and voltages on y axis and the data as the heatmap
# plt.figure(figsize=(30, 6))
# sns.heatmap(data, cmap='gist_rainbow')
# plt.title('Heatmap of Raman Data')
# plt.ylabel('Voltage (V)')
# plt.xlabel('Wavenumber (cm^-1)')
# plt.xticks(range(len(x_axis)), x_axis, rotation=90)  # Rotate x-axis labels by 90 degrees
# plt.tight_layout()  # Adjust spacing to prevent overlapping labels
# plt.savefig(f'heatmap_{sheet}.png', dpi=1200)
# print('data shape', data.shape)

# # Do the same thing as above but with imshow
# plt.figure(figsize=(10, 6))
# plt.imshow(data, cmap='viridis', interpolation='spline16')
# plt.colorbar()
# plt.title('Heatmap of Raman Data')
# plt.ylabel('Voltage (V)')
# plt.xlabel('Wavenumber (cm^-1)')
# # plt.xticks(range(len(x_axis)), x_axis, rotation=90)  # Rotate x-axis labels by 90 degrees
# plt.savefig('heatmap_matplotlib.png', dpi=300)

# Create a contour plot
plt.figure(figsize=(10, 6))
plt.contourf(x_axis,y_axis, data, cmap='viridis')
plt.colorbar()
plt.title('Contour Plot of Raman Data')
plt.ylabel('Voltage (V)')
plt.xlabel('Wavenumber (cm^-1)')
plt.savefig(f'contour_{sheet}.png', dpi=300)
