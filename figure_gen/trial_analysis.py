import seaborn as sns
import matplotlib.pyplot as plt
import pandas
# csv = pandas.read_csv('../.csv')
# filename = 
import os
import numpy as np
from numpy.polynomial import polynomial as P
import csv

directory = os.fsencode('../plots')
jerk_norms = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".csv"):
        print(filename)
        csv = pandas.read_csv('../plots/'+filename)
        
        # df = sns.load_dataset("titanic")
        fig = plt.gcf()
        fig.clear()
        fig.set_size_inches(7, 10)
        plt.subplot(2,2,1)
        print(csv.keys())
        plt.plot(csv['t'], csv['x_pos'])
        coeff_px = np.polyfit(csv['t'], csv['x_pos'], 99)
        val_px   = np.polyval(coeff_px, csv['t'])
        print(csv['t'].count())
        plt.plot(csv['t'], val_px)

        coeff_vx = np.polyder(coeff_px)
        val_vx   = np.polyval(coeff_vx, csv['t'])
        plt.subplot(2,2,2)
        plt.plot(csv['t'], csv['x_vel'])
        plt.plot(csv['t'], val_vx)
        # these verified taking derivative after polyfit is ok!!

        plt.subplot(2,2,3)
        coeff_vx_data = np.polyfit(csv['t'], csv['x_vel'], 99)
        coeff_vy_data = np.polyfit(csv['t'], csv['y_vel'], 99)
        coeff_vz_data = np.polyfit(csv['t'], csv['z_vel'], 99)
        coeff_jx = np.polyder(coeff_vx_data, m = 2)
        coeff_jy = np.polyder(coeff_vy_data, m = 2)
        coeff_jz = np.polyder(coeff_vz_data, m = 2)
        val_jx   = np.polyval(coeff_jx, csv['t'])
        val_jy   = np.polyval(coeff_jy, csv['t'])
        val_jz   = np.polyval(coeff_jz, csv['t'])
        coeff_jx2   = P.polypow(coeff_jx[::-1], 2)
        
        plt.subplot(2,2,4)
        val_jx2   = P.polyval(csv['t'], coeff_jx2)
        plt.plot(csv['t'], val_jx)
        plt.plot(csv['t'], val_jx2) #squaring 100 degree polynomial results in big problems things blow up
        plt.plot(csv['t'], (val_jx**2))
        plt.ylim(-1,6)
        plt.legend(['jerk', 'jerk squared polyval', 'jerk squared from jerk**2'])
        jerk_norm_avg_time = np.sqrt(np.trapz(val_jx**2+val_jy**2+val_jz**2, csv['t']))/csv.iloc[-1]['t']
        jerk_norms.append(jerk_norm_avg_time)
        print("2 norm of jerk", jerk_norm_avg_time)

        # coeff_jx

        
        # plt.plot(csv['t'], csv['y_pos'])
        # plt.plot(csv['t'], csv['z_pos'])
        # plt.show()
        # plt.subplot(2,2,2)
        # sns.boxplot(x=csv['u_norm_over_hover'],showfliers = True)
        # print(csv[["u_smooth"]].describe())
        # plt.xlim(0, 10)
        # # plt.xlim(0, 2)
        # plt.subplot(2,2,3)
        # sns.boxplot(x=csv['traj_cost'],showfliers = True)
        # plt.xlim(-0.5,1.5)
        # plt.subplot(2,2,4)
        # sns.boxplot(x=csv['solve_time'],showfliers = True)
        # plt.xlim(0,400)
        # plt.suptitle(filename)
        # groupby = csv.groupby('solver_status').count()
        # print(groupby)

        # plt.show()
        # fig.savefig(filename+".png", dpi=100)

