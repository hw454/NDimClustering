# Crop all of the matrices in mat_list
#
# :description: include only the first n_rows
#   and columns between n_col0 and n_col1,
#   ensure that out_pheno and exp_pheno are included.
#
#   1. Get the columns of each matrix between n_col0 and n_col1.
#     Crop to n_rows at this step at the same time
#   2. Check if the exposure and outcome columns are in the n_col range.
#   3. If not then:
#      a. Get the matrix with just the exposure and outcome columns
#      using `crop_mat_colnames`. Crop to n_rows at this step at the same time
#   4. Merge the two sets to create a matrix with columns in the n_col
#   range and the exposure and outcome columns.
#
# :param: mat_list The list of matrices that need to be cropped.
# :param: trait_df The dataframe of the traits.
# :param: out_pheno The outcome phenotype label
# :param: exp_pheno The ecposure phenotype label.
# :param: n_rows The number of rows to crop by.
# :param: n_col0 The first column to include
# :param: n_col1 The last column to include.
#
# :return: List of cropped matrices.
#
# @family preconditioning_functions
import pandas as pd
def crop_matrices (mat_dict, trait_df,
  out_pheno, exp_pheno, n_rows, n_col0, n_col1
) :
  mat_names = mat_dict.keys
  mat_nrows = mat_dict[mat_names[0]]

  # If the data doesn't have enough rows to match the test request return all available rows.
  if (n_rows > mat_nrows):
    n_rows <- mat_nrows
  
  mat_dict_out = {}
  for key in mat_names:
    mat_dict_out[key] = mat_dict[key].iloc[0:n_rows, n_col0:n_col1]

  out_check =  out_pheno in trait_df.phenotype[n_col0:n_col1]
  exp_check = exp_pheno in trait_df.phenotype[n_col0:n_col1]

  # Check if the desired outcome and exposure are in the test columns
  if (not out_check or not exp_check) :

    mat_out_exp = {}
    # For each matrix crop the rows and get the exposure and outcome columns
    for key in mat_names:
      mat_out_exp[key] = mat_dict[key].iloc[0:n_rows]
      mat_out_exp[key] = mat_dict[key][[out_pheno, exp_pheno]]
      # Combine matrices and remove duplicates
      mat_dict_out[key] = pd.concat([mat_out[key], mat_out_exp[key]], axis =1).drop_duplicates()
    trait_list = mat_dict_out[mat_names[0]].columns()
    trait_df.phenotype = trait_list
    # Collect the matrices into one object
    data_matrices = mat_dict_out,+{"trait_info" : trait_df}
  else :
    trait_list = mat_dict_out[mat_names[0]].columns()
    data_matrices = mat_dict_out + {"trait_info" : trait_list}
  return(data_matrices)
