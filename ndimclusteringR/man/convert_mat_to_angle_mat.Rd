% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_mat_to_angle_mat.R
\name{convert_mat_to_angle_mat}
\alias{convert_mat_to_angle_mat}
\title{Convert a matrix of spatial data into a matrix of angles}
\usage{
convert_mat_to_angle_mat(mat)
}
\arguments{
\item{mat}{the matrix of spatial data, each row corresponds to a point.}
}
\value{
ang_mat
}
\description{
The angles correspond to the angles to the unit vector
  on each axis.
  The angles are calculated for each point (row) in mat
  using [convert_point_to_angles_all] which uses
  [convert_point_to_angle] for each column.
  The angles for each point are then row bound using "rbind".
  into "ang_mat".
}
\seealso{
Other clustering_components: 
\code{\link{calc_clust_cent}()},
\code{\link{calc_clust_prob}()},
\code{\link{calc_member_dist_cent}()},
\code{\link{check_clust_cent}()},
\code{\link{get_aic}()},
\code{\link{km_nan}()},
\code{\link{make_clust_col_name}()},
\code{\link{make_clust_col}()},
\code{\link{make_rand_cent}()}

Other preconditioning_functions: 
\code{\link{convert_point_to_angle}()},
\code{\link{crop_mat_colnames}()},
\code{\link{crop_mat_colnums}()},
\code{\link{crop_matrices}()},
\code{\link{reformat_data}()},
\code{\link{rescale_by_end_col}()}
}
\concept{clustering_components}
\concept{preconditioning_functions}
