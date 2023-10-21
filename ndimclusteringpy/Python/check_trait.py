import pandas as pd
""" Check is a trait is valid

 :param a: The string label for the trait.
 :param pheno_list: The list of all traits.
 :param data_df: The dataframe containing the data of interest.
   With columns labelled by trait labels.
 :param na_percent: percentage of a column which is required to be NOT NaN.

  .. code-block:: python
   
     if trait passes [na_col_check]:
       trait_single_df = data.frame(label = a,
                              axes_ind = position in pheno_list)
     else: trait_single_df = data.frame()

 :return: trait_single_df
 :rtype: Dataframe"""
def check_trait(a, pheno_list, data_df, na_percent):
  # Add the trait to the trait dataframe
  allna = check_col_na(data_df[a], na_percent)
  if (not allna) :
    a_ind = np.where(pheno_list == a)[0]
    trait_single_df = pd.DataFrame(data = {
                                   "label" : a,
                                   "axes_ind" : a_ind})
  else:
    trait_single_df = pd.DataFrame(columns = [
                                    "label",
                                    "axes_ind"])
  return(trait_single_df)