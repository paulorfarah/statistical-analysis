# Required imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from roses.statistical_test.kruskal_wallis import kruskal_wallis
from roses.effect_size import vargha_delaney

# rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# For a beautiful plot
plt.style.use('ggplot')

# If you want let the plot beautiful
import seaborn as sns
sns.set_style("whitegrid")
sns.set(palette="pastel")

# Improving the readability in the plots
FONT_SIZE_PLOTS = 16
plt.rcParams.update({
    'font.size': FONT_SIZE_PLOTS,
    'xtick.labelsize': FONT_SIZE_PLOTS,
    'ytick.labelsize': FONT_SIZE_PLOTS,
    'legend.fontsize': FONT_SIZE_PLOTS,
    'axes.titlesize' : FONT_SIZE_PLOTS,
    'axes.labelsize' : FONT_SIZE_PLOTS
})

def print_mean(df, metric):
    """
    This function group the data and presents for each algorithm the mean, std, max, and min values.
    Besides that, it returns the best algorithm based on the maximum mean value.
    """
    mean = df.groupby(['Application'], as_index=False).agg({metric: ['mean', 'std', 'max', 'min']})
    mean.columns = ['name', 'mean', 'std', 'max', 'min']
    
    # Round values (to be used in the article)
    mean = mean.round({'mean': 4, 'std': 3, 'max': 4, 'min': 4})
    mean = mean.infer_objects()

    bestalg = mean.loc[mean['mean'].idxmax()]

    return mean, bestalg['name']

if __name__ == '__main__':
    DATASETS = ['results-methods-applications.csv']
    metrics = ['F1', 'Acc', 'Sen', 'AUC']
    models = ['Str', 'StrEvo', 'Evo']
    #models = ['RUS','ROS','ENN','ADA','SMOTE','NONE']

    for dataset in DATASETS:
        for metric in metrics:
            for model in models:
                df = pd.read_csv('resources/' + dataset, sep=";")
   
                #df.drop(indexNames , inplace=True)
                fig, ax = plt.subplots(figsize=(10, 6))
               
                # for KRUSKAL WALLIS ANALYSIS
                k = kruskal_wallis(df, metric, 'Application', sort=False)
                kruskal_result, posthoc = k.apply(ax)


                
                plt.tight_layout()
                #kruskal_results
                # ax.set_ylabel(model + " - " + metric )
                ax.set_ylabel(metric)
                plt.savefig("results/results-methods-applications-" + metric + "-" + model + ".pdf", bbox_inches='tight')
                plt.cla()
                plt.close(fig)
                plt.close()

                mean, best = print_mean(df, metric)

                print("Best config:", best, "\n\nMeans:")
                try:
                    #print(posthoc[0])
                    print("***********")
                    print(posthoc[0])
                    print("***********")
                    print(posthoc[1])
                    print("***********")
                    # Get the posthoc and prepare the comparison against the best algorithm
                    df_eff = vargha_delaney.reduce(posthoc[1], best)

                    # Define symbols for each effect size magnitude
                    mean['eff_symbol'] = " "

                    def get_eff_symbol(x, best, df_eff):
                        if x['name'] == best:
                            return "$\\bigstar$"
                        elif len(df_eff.loc[df_eff.compared_with == x['name'], 'effect_size_symbol'].values) > 0:
                            return df_eff.loc[df_eff.compared_with == x['name'], 'effect_size_symbol'].values[0]
                        else:
                            return df_eff.loc[df_eff.base == x['name'], 'effect_size_symbol'].values[0]

                    mean['eff_symbol'] = mean.apply(lambda x: get_eff_symbol(x, best, df_eff), axis=1)

                    # Concat the values to a unique columns
                    mean['avg_std_effect'] = mean.apply(lambda row: f"{row['mean']:.4f} $\\pm$ {row['std']:.4f} {row['eff_symbol']}", axis=1)

                    # Select the main information
                    mean = mean[['name', 'avg_std_effect']]

                    print(mean)

                    mean_trans = mean.copy()
                    mean_trans.index = mean['name']
                    mean_trans = mean_trans.transpose()

                    temp_x = mean_trans.to_string(index=False, index_names=False).split("\n")[1:]
                    print(temp_x[0]) # Column names

                    # Just a beautiful print to use with LaTeX table :}
                    temp_split = temp_x[1].split()
                    cont = 1
                    for x in temp_split:
                        print(f"{x} ", end="")
                        if (cont % 4 == 0 and cont != len(temp_split)):
                            print("& ", end="")
                        cont += 1
                    print("\n\n")

                    ro.r.assign('posthoc', posthoc[0])
                    ro.r('posthoc_table <- t(as.matrix(posthoc$p.value))')
                    ro.r('df_posthoc <- as.data.frame(t(posthoc_table))')

                    with localconverter(ro.default_converter + pandas2ri.converter):
                        posthoc_df = ro.conversion.rpy2py(ro.r('df_posthoc'))

                    def get_config_latex(row, best, posthoc_df):
                        """
                        Latex commands used:
                        - Best algorithm: \newcommand{\cellgray}[1]{\cellcolor{gray!30}{#1}}
                        - Equivalent to the best one: \newcommand{\cellbold}[1]{\cellcolor{gray!30}{\textbf{#1}}}
                        """
                        current_name =  row['name']
                        is_equivalent = False
                        #print(current_name)
                        #print(best)
                        #print(posthoc_df.head(15))
                        
                        if best in posthoc_df.columns and current_name in posthoc_df.index and not np.isnan(posthoc_df.loc[current_name][best]):
                            # They are equivalent
                            is_equivalent = posthoc_df.loc[current_name][best] >= 0.05                    
                        elif current_name in posthoc_df.columns and best in posthoc_df.index and not np.isnan(posthoc_df.loc[best][current_name]):
                            # They are equivalent
                            is_equivalent = posthoc_df.loc[best][current_name] >= 0.05             
                            
                        if is_equivalent:
                            return f"\\cellgray{{{row['avg_std_effect']}}}"
                        elif row['name'] == best:
                            return f"\\cellbold{{{row['avg_std_effect']}}}"
                        else:
                            return row['avg_std_effect'] 
                                    
                    mean['latex_format'] = mean.apply(lambda row: get_config_latex(row, best, posthoc_df), axis=1)

                        # Select the main information
                    
                    mean = mean[['name', 'latex_format']]

                    mean_trans = mean.copy()
                    mean_trans.index = mean['name']
                    mean_trans = mean_trans.transpose()

                    temp_x = mean_trans.to_string(index=False, index_names=False).split("\n")[1:]
                    #print(temp_x[0]) # Column names
                    #print(np.matrix(temp_x))
                    
                    # Just a beautiful print to use with LaTeX table :}
                    temp_split = temp_x[1].split()
                    cont = 1
                    for x in temp_split:
                        print(f"{x} ", end="")
                        if (cont % 4 == 0 and cont != len(temp_split)):
                            print("& ", end="")
                        cont += 1
                    print("\n\n")

                except:
                    print("Sem diferenca estatistica metrica " + metric)

            
