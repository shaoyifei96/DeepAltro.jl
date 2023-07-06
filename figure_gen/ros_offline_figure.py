import seaborn as sns
import matplotlib.pyplot as plt
import pandas
# csv = pandas.read_csv('../.csv')
# filename = 
import os

directory = os.fsencode('./')
    
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".csv") and filename.startswith("deleting"):
        print(filename)
        csv = pandas.read_csv(filename)
        print("iterations = ",csv.loc[:, 'solver_iter'].mean())
        # df = sns.load_dataset("titanic")
        fig = plt.gcf()
        fig.clear()
        fig.set_size_inches(7, 10)
        plt.subplot(2,2,1)
        sns.violinplot(x=csv['u_norm_over_hover'])
        plt.subplot(2,2,2)
        sns.boxplot(x=csv['u_norm_over_hover'],showfliers = False)
        plt.xlim(0, 2)
        plt.subplot(2,2,3)
        sns.violinplot(x=csv['solve_time'])
        plt.subplot(2,2,4)
        sns.boxplot(x=csv['solve_time'],showfliers = False)
        plt.xlim(15,80)
        plt.suptitle(filename)
        groupby = csv.groupby('solver_status').count()
        print(groupby)

        # plt.show()
        fig.savefig(filename+".png", dpi=100)

