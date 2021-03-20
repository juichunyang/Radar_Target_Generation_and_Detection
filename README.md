# Radar Target Generation and Detection

This is the third project of Sensor Fusion Nanodegree of Udacity. I calculated the Range-Doppler Map (RDM) to ascertain the position and velocity of the target. The following picture is the layout of this project.

<img src="https://github.com/CuteJui/Radar_Target_Generation_and_Detection/blob/main/readme_resource/layout.png" width="800">

## Radar Specifications
- Frequency of operation : 77 GHz
```
fc = 77e9 
```
- Max Range : 200 m
```
maxRange = 200
```
- Range Resolution : 1 m
```
rangeResolution = 1
```
- Max Velocity : 100 m/s 
```
maxV = 100
```

## CFAR Parameters
- Number of Training Cells in both the dimensions
```
Tr = 12
Td = 10
```
- Number of Guard Cells in both the dimensions
```
Gr = 5
Gd = 5
```
- Offset the threshold by SNR value in dB
```
offset = 6
```
