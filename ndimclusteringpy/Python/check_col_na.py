""" Check the percentage of NaNs in column

 :param percent: The percentage of the column which must be not NaN
 :param b_col: The column vector from the matrix

 :return: 0 or 1 - 1 if `b_col` has an acceptable number of non_NaN, else 0.
 :rtype: Bool"""
def check_col_na(b_col, percent=0.95):
    n_accept = len(b_col) * (1 - percent)
    narows = b_col.isna().sum()
    if len(narows) > n_accept:
        return True
    else:
        return False