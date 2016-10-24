/*** quality_measures.cpp
 * version 0.1
 * (c) 2016 Gregory Breard, University of Rhode Island
 *
 * This file contains a set of functions used for 
 * evaluating the quality of self-organizing maps (SOMs).
 *** License
 * This program is free software; you can redistribute it 
 * and/or modify it under the terms of the GNU General 
 * Public License as published by the Free Software 
 * Foundation.
 *
 * This program is distributed in the hope that it will 
 * be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * A copy of the GNU General Public License is available 
 * at: http://www.r-project.org/Licenses/
 ***/

#include <Rcpp.h>

using namespace Rcpp;

// Reference:
// T. Kohonen, Self-organizing maps, Berlin: Springer, 
// 2001.
// [[Rcpp::export(name = "get.quant.err")]]
List GetQuantizationError(NumericMatrix dist_cross) {
  // Initialize
  int n = dist_cross.nrow();
  double total_dist = 0;
  
  // Check each data point
  for (int i = 0; i < n; i++) {
    NumericVector between = dist_cross.row(i);
    double bmu_dist = min(between);
    total_dist += bmu_dist;	
  } // end for (i)
  
  // Calculate the error
  double err = total_dist / n;

  // Return list
  List out = List::create(Named("val") = err);

  return out;
} // end GetQuantizationError

// Reference:
// G. Polzlbauer, Survey and comparison of quality 
// measures for self-organizing maps, in Proc. 5th 
// Workshop Data Analysis, pg 67–82, 2004.
// [[Rcpp::export(name = "get.top.err")]]
List GetTopographicError(NumericMatrix dist_cross, 
                         int xdim) {
  // Initialize
  int n = dist_cross.nrow();
  int errors = 0;
  
  // Check each data point
  for (int i = 0; i < n; i++) {
    // Intialize variables
    NumericVector between = dist_cross.row(i);
    
    // Initialize the unsorted index vectors
    std::vector<int> idx(between.size());
    std::iota(idx.begin(), idx.end(), 0);
    
    // Sort the indices by the distance
    std::sort(idx.begin(), idx.end(),
              [between](double i1, double i2) {
                    return between[i1] < between[i2];
                });
    
    // Get the best (and second best) matching units
    int bmu_index = idx[0];
    int sbmu_index = idx[1];
    
    // Considering the neighborhood:
    // n-xdim-1 n-xdim n-xdim+1
    //   n-1       n      n+1
    // n+xdim+1 n+xdim n+xdim+1
    
    // Find index difference
    int dif = abs(bmu_index - sbmu_index);
    
    // Check for error
    if (!(dif == 1 || dif == xdim - 1 
                          || dif == xdim 
                              || dif == xdim + 1))
        errors++;	
  } // end for (i)
  
  // Calculate the error
  double err = (double)errors / n;
  
  // Return list
  List out = List::create(Named("val") = err);

  return out;
} // end GetTopographicError

// Reference:
// T.  Villmann, R. Der, M. Herrmann, and T. Martinetz, 
// Topology preservation in self-organizing feature maps:
// exact definition and measurement, IEEE Trans. Neural 
// Netw., vol. 8 no. 2, pg 256 - 266, 1997.
// [[Rcpp::export(name = "get.top.func")]]
List GetTopographicFunction(NumericMatrix dist_cross, 
                            int xdim) {
  // Initialize
  int n = dist_cross.nrow();
  int m = dist_cross.ncol();
  
  // Initialize the connectivity and Delaunay 
  // Triangulation matrices
  NumericMatrix C(m, m);
  
  // Build the connectivity matrix
  for (int i = 0; i < n; i++) {
    // Intialize variables
    NumericVector between = dist_cross.row(i);
    
    // Initialize the unsorted index vectors
    std::vector<int> idx(between.size());
    std::iota(idx.begin(), idx.end(), 0);
    
    // Sort the indices by the distance
    std::sort(idx.begin(), idx.end(),
              [between](double i1, double i2) {
                      return between[i1] < between[i2];
                });
    
    // Get the best (and second best) matching units
    int bmu_index = idx[0];
    int sbmu_index = idx[1];
    
    // Add an edge between the best and second 
    // best matching units
    C(bmu_index, sbmu_index) = 1;
    C(sbmu_index, bmu_index) = 1;	
  } // end for (i)
  
  // Build the Delaunay Triangulation matrix (shortest
  // paths) using Floyd–Warshall algorithm
  //   initialize paths
  NumericMatrix Dm(clone(C));
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      if (i == j)
        Dm(i, j) = 0;
      else if (Dm(i, j) != 1)
        Dm(i, j) = std::numeric_limits<double>::infinity();
  //   find shortest paths
  for (int k = 0; k < m; k++)
    for (int i = 0; i < m; i++)
      for (int j = 0; j < m; j++)
        if (Dm(i, k) + Dm(k, j) < Dm(i, j))
          Dm(i, j) = Dm(i, k) + Dm(k, j);
  
  // Initialize function results
  NumericVector ks = NumericVector(2 * m - 1);
  NumericVector phi = NumericVector(2 * m - 1);
  
  // Check that we have a valid triangulation
  if (max(Dm) == std::numeric_limits<double>::infinity()) {
    // Can't calculate
    for (int i = 0; i < 2 * m - 1; i++) {
      ks(i) = i - m + 1;
      phi(i) = nan("");
    } // end for (i)
  } else {
    // Calculates all function values
    for (int i = 0; i < 2 * m - 1; i++) {
      int k = i - m + 1;
      double p = 0.0;
      
      // Calculate phi(k)
      for (int j = 0; j < m; j++) {
        for (int l = 0; l < m; l++) {
          int f = 0;
          NumericVector i_idx(2);
          NumericVector j_idx(2);
          i_idx(0) = j % xdim;
          i_idx(1) = floor(j / xdim);
          j_idx(0) = l % xdim;
          j_idx(1) = floor(l / xdim);
        
          // Calculate f(k)
		  double dist_Dm = Dm(j, l);
          if (k > 0) {
            double dist = max(abs(i_idx - j_idx));
            if (dist > k && dist_Dm == 1)
              f++;
          } else if (k < 0) {
            double dist = sum(abs(i_idx - j_idx));
            if (dist == 1 && dist_Dm > abs(k))
              f++;
          } // end if
          
          p = p + f;
        } // end for (l)
      } // end for (j)
  
      ks(i) = k;
      phi(i) = p / m;
    } // end for (i)
    
    // Set phi(0)
    phi(m) = phi(m + 1) + phi(m - 1); 
  } // end if
  
  // Return list
  List out = List::create(Named("k") = ks,
                          Named("phi") = phi);
  
  return out;
} // end GetTopographicFunction

// Reference:
// J. Venna and S. Kaski, "Neighborhood preservation in 
// nonlinear projection methods: An experimental study", 
// Lecture Notes in Comput. Sci., vol. 2130, pg 485-491, 
// 2001.
// [[Rcpp::export(name = "get.hood.pres")]]
List GetNeighborhoodPreservation(NumericMatrix dist_data, 
                                 NumericMatrix dist_proj, 
                                 int k) {
  // Initialize
  int n = dist_data.nrow();
  double M_1 = 0.0;
  double M_2 = 0.0;
  
  // Check each data point
  for (int i = 0; i < n; i++) {
    // Get the distances for x_i
    NumericVector dist = dist_data.row(i);
    NumericVector pdist = dist_proj.row(i);
    
    // Initialize the unsorted index vectors
    std::vector<int> idx(dist.size());
    std::vector<int> pidx(pdist.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::iota(pidx.begin(), pidx.end(), 0);
    
    // Sort the indices by the distance
    std::sort(idx.begin(), idx.end(),
                [dist](double i1, double i2) {
                      return dist[i1] < dist[i2];
                  });
    std::sort(pidx.begin(), pidx.end(),
                [pdist](double i1, double i2) {
                      return pdist[i1] < pdist[i2];
                  });
    
    // Get all x_j in (and not in)  C_k(x_i) 
    std::vector<int> Ck(idx.begin(), idx.begin() + k);
    std::vector<int> not_Ck(idx.begin() + k, idx.end());
    
    // Get all x_j in (and not in) C^_k(x_i)
    std::vector<int> hat_Ck(pidx.begin(), 
                            pidx.begin() + k);
    std::vector<int> not_hat_Ck(pidx.begin() + k, 
                                pidx.end());
    
    // Get U_k(x_i), e.g. the intersection of x_j 
    // not in C_k(x_i) and x_j in C^_k(x_i)
    std::vector<int> Uk(k);
    std::sort(not_Ck.begin(), not_Ck.end());
    std::sort(hat_Ck.begin(), hat_Ck.end());
    std::vector<int>::iterator it 
        = std::set_intersection(not_Ck.begin(), 
                                not_Ck.end(), 
                                hat_Ck.begin(), 
                                hat_Ck.end(), 
                                Uk.begin());
    Uk.resize(it - Uk.begin());
    
    // Get V_k(x_i), e.g. the intersection of x_j in
    // C_k(x_i) and x_j not in C^_k(x_i)
    std::vector<int> Vk(k);
    std::sort(Ck.begin(), Ck.end());
    std::sort(not_hat_Ck.begin(), not_hat_Ck.end());
    it = std::set_intersection(Ck.begin(), Ck.end(), 
                               not_hat_Ck.begin(), 
                               not_hat_Ck.end(), 
                               Vk.begin());
    Vk.resize(it - Vk.begin());
    
    // Calculate the inner sums
    for (int j = 0; j < std::max(Uk.size(), 
                                 Vk.size()); j++) {
      if (j < Uk.size()) {
          int x_j = Uk[j];
          if (x_j != i) {
            int r = find(idx.begin(), idx.end(), x_j) 
                                - idx.begin() + 1;
            M_1 = M_1 + r - k;
          } // end if (x_j)
      } // end if (j)
      if (j < Vk.size()) {
        int x_j = Vk[j];
        if (x_j != i) {
          int r_hat = find(pidx.begin(), pidx.end(), x_j)
                                - pidx.begin() + 1;
          M_2 = M_2 + r_hat - k;
        } // end if (x_j)
      } // end if (j)
    } // end for (j)
  } // end for (i)
  
  // Convert the sum
  M_1 = 1 - (2 * M_1 / (n * k * (2 * n - 3 * k - 1)));
  M_2 = 1 - (2 * M_2 / (n * k * (2 * n - 3 * k - 1)));
  
  // Return list
  List out = List::create(Named("k") = k,
               Named("trustworthiness") = M_1,
               Named("neighborhood.preservation") = M_2);
  
  return out;
} // end GetNeighborhoodPreservation

// Reference:
// J. Hirschberg and A. Rosenberg, V-Measure: A 
// conditional entropy-based external cluster evaluation,
// Columbia University Academic Commons, 2007, 
// http://hdl.handle.net/10022/AC:P:21139
// [[Rcpp::export(name = "get.v.measure")]]
List GetVMeasure(IntegerVector labels, 
                 IntegerVector clusters, 
                 double beta = 1.0) {
  // Check the sizes
  if (labels.size() != clusters.size())
    stop("get.v.measure: vectors sizes don't match.");
  
  // Get the level sizes
  int N = labels.size();
  int n = sort_unique(labels).size();
  int m = sort_unique(clusters).size(); 
  
  // Generate the contingency table
  NumericMatrix A(n, m);
  for (int i = 0; i < N; i++) {
    int l = labels[i] - 1;
    int c = clusters[i] - 1;
    A(l, c) = A(l, c) + 1;
  }
  //   convert to probabilities for entropy (H)
  A = A / N;

  // Initialize values
  double H_CK = 0.0;
  double H_C = 0.0;
  double H_KC = 0.0;
  double H_K = 0.0;
  double homo = 0.0;
  double comp = 0.0;
  
  // Calculate H(C|K)
  for (int k = 0; k < m; k++)
    for (int c = 0; c < n; c++)
      if (A(c, k) > 0)
        H_CK = H_CK + A(c, k) * (log(A(c, k)) 
                                   - log(sum(A(_, k))));
  H_CK = - H_CK;
  
  // Calculate H(C)
  for (int c = 0; c < n; c++)
    H_C = H_C + sum(A(c, _)) * log(sum(A(c, _)));
  H_C = - H_C;
  if (std::isnan(H_C))
    H_C = 0;
  
  // Calculate H(K|C)
  for (int c = 0; c < n; c++)
    for (int k = 0; k < m; k++)
      if (A(c, k) > 0)
        H_KC = H_KC + A(c, k) * (log(A(c, k)) 
                                   - log(sum(A(c, _))));
  H_KC = - H_KC;
  
  // Calculate H(K)
  for (int k = 0; k < m; k++)
    H_K = H_K + sum(A(_, k)) * log(sum(A(_, k)));
  H_K = - H_K;
  if (std::isnan(H_K)) H_K = 0;
  
  // Calculate homogeneity
  if (H_C == 0) homo = 1;
  else homo = 1 - H_CK / H_C;
  
  // Calculate completeness
  if (H_K == 0) comp = 1;
  else comp = 1 - H_KC / H_K;
 
  // Calculate the weighted harmonic mean of 
  // homogeneity and completeness
  double v = ((1 + beta) * homo * comp) / 
                              ((beta * homo) + comp);
  
  // Return list
  List out = List::create(Named("beta") = beta,
                          Named("homogeneity") = homo,
                          Named("completeness") = comp,
                          Named("v.measure") = v);
  
  return out;
} // end GetVMeasure
