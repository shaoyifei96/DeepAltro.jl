import seaborn as sns
import matplotlib.pyplot as plt
import pandas
# csv = pandas.read_csv('../.csv')
# filename = 
import os

directory = os.fsencode('./')
    
two_file_comp = True
pandas.set_option('display.max_rows', 10)
if two_file_comp:
    file1 = "tilos_no_cons_hover.csv"
    file2 = "tilos_with_cons_hover.csv"
    csv1 = pandas.read_csv(file1)
    csv2 = pandas.read_csv(file2)
    #show a few rows of each
    print(csv1.head())
    print(csv2.head())
    print(csv1.describe())
    print(csv2.describe())

    #subtract u_smooth column from both cvs and print the difference
    res1 =  csv2['u_smooth'] - csv1['u_smooth']
    print("most reduced u smoothness in trial", res1.idxmin())
    print("u-smooth", res1)
    print(res1.describe())
    #do the same for u_norm_over_hover column
    res2 =  csv2['u_norm_over_hover'] - csv1['u_norm_over_hover']
    res3 = csv2['solve_time'] - csv1['solve_time']
    print("most reduced effort in trial", res2.idxmin())
    print("u_norm_over_hover", res2)
    print(res2.describe())
    print(res3.describe())
    plt.subplot(2,2,1)
    sns.boxplot(res1)
    plt.subplot(2,2,2)
    sns.boxplot(res2)
    plt.show()



else:

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and filename.startswith("tilos"):
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

