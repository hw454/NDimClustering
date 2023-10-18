""" Add "np" columns to a dataframe.

  Representing the principal components a column for the PC score
  is added for each PC.
  The label for each column is given by the :func:"make_p_cols".

 :param df: the dataframe to add columns to
 :param np: the number of columns to add.

 :return: df
 :rtype dataframe
"""
def add_np_cols (df, np):
  # Create dataframe with columns with labels Pi
  # for i in 1 to np to the dataframe.
  p_cols = [make_p_col(i) for i in range(1,np)]
  np_df = pd.concat(p_cols, axis=1)
  full_df = pd.concat([df, np_df], axis=1)
  return(full_df)