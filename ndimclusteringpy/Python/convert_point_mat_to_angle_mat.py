import numpy as np
import pandas as pd
""" Convert a matrix of spatial data into a matrix of angles

 :param df: the dataframe matrix of spatial data, each row corresponds to a point.

 The angles correspond to the angles to the unit vector
   on each axis.
   The angles are calculated for each point (row) in mat
   using :func:`convert_point_to_angles_all` which uses
   :func:`convert_point_to_angle` for each column.
   The angles for each point are then row bound using pandas concat.
   into "ang_mat".

 :return: ang_mat"""
def mat_to_angle_mat (df) :
  nr, nc = df.shape
  unit_mat = np.identiy(nc)
  p_list = [convert_point_to_angles_all(j,
                  mat = df,
                  unit_mat = unit_mat) 
                  for j in range(1,nr)]
  ang_mat = pd.concat(p_list)
  ang_mat.index = df.index
  ang_mat.columns = df.columns
  return(ang_mat)