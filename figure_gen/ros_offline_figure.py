import seaborn as sns
import matplotlib.pyplot as plt
import pandas
# csv = pandas.read_csv('../.csv')
# filename = 
import os

directory = os.fsencode('./')
    
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    print(filename)
    if filename.endswith(".csv") and filename.startswith("new_test_cons"):
        print(filename)
        csv = pandas.read_csv(filename)
        print("iterations = ",csv.loc[:, 'solver_iter'].mean())
        # df = sns.load_dataset("titanic")
        fig = plt.gcf()
        fig.clear()
        fig.set_size_inches(7, 10)
        plt.subplot(2,2,1)
        sns.boxplot(x=csv['u_smooth'],showfliers = True)
        csv[["u_smooth"]].describe()
        plt.xlim(0, 13)
        plt.subplot(2,2,2)
        sns.boxplot(x=csv['u_norm_over_hover'],showfliers = True)
        print(csv[["u_smooth"]].describe())
        plt.xlim(0, 10)
        # plt.xlim(0, 2)
        plt.subplot(2,2,3)
        sns.boxplot(x=csv['traj_cost'],showfliers = True)
        plt.xlim(-0.5,1.5)
        plt.subplot(2,2,4)
        sns.boxplot(x=csv['solve_time'],showfliers = True)
        plt.xlim(0,400)
        plt.suptitle(filename)
        groupby = csv.groupby('solver_status').count()
        print(groupby)

        # plt.show()
        fig.savefig(filename+".png", dpi=100)

