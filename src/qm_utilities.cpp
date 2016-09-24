#ifndef __QM_UTILITIES__
#define __QM_UTILITIES__

#include <Rcpp.h>
using namespace Rcpp;

/**
* Converts a DataFrame into a NumericMatrix.
* @param x The data frame to be converted.
* @return The matrix.
*/
inline NumericMatrix df_to_mat(const DataFrame x) {
  // Initialize matrix
  NumericMatrix y(x.nrows(), x.size());
  
  // Copy df rows to matrix
  for (int i = 0; i < x.size(); i++)
    y(_, i) = NumericVector(x[i]);
  
  return y;
} // end df_to_mat

/**
 * Gets the Euclidean distance between two NumericVectors.
 * @param x The first vector.
 * @param y The second vector.
 * @return The distance.
 */
inline double calc_dist(const NumericVector x, const NumericVector y) {
  NumericVector inside = pow(x - y, 2);
  return sqrt(sum(inside));
} // end dist

/**
 * Gets the Euclidean distances between all vectors in a NumericMatrix.
 * @param x The matrix.
 * @return The distance matrix.
 */
inline NumericMatrix calc_dist_mat (const NumericMatrix x){
  // Initialize distance matrix
  NumericMatrix dist_mat(x.nrow(), x.nrow());
  
  // Find the distance between all pairs
  for (int i = 0; i < x.nrow() - 1; i++){
    NumericVector v1 = x.row(i);
    
    for (int j = i + 1; j < x.nrow(); j ++){
      double dist = calc_dist(v1, x.row(j));
      dist_mat(j,i) = dist;
      dist_mat(i,j) = dist;
    } // end for (j)
  } // end for(i)
  
  return dist_mat;
} // end calc_dist_mat

/**
 * Gets the max distance between two NumericVectors.
 * @param x The first vector.
 * @param y The second vector.
 * @return The distance.
 */
inline double calc_dist_max(const NumericVector x, const NumericVector y) {
  return max(x - y);
} // end dist


#endif // __QM_UTILITIES__