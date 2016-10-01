#include <cmath>
#include "qm_utilities.cpp"

using namespace Rcpp;

// Reference:
// T. Kohonen, Self-organizing maps, Berlin: Springer, 2001.
// [[Rcpp::export(name = "get.quant.err")]]
List GetQuantizationError(List map) {
  // Create matrices for the data and neurons
  NumericMatrix data = df_to_mat(DataFrame(map["data"]));
  NumericMatrix neurons = df_to_mat(DataFrame(map["neurons"]));
  
  // Initialize distance
  double total_dist = 0;
  
  // Check each data point
  for (int i = 0; i < data.nrow(); i++) {
    // Intialize variables
    double bmu_dist = std::numeric_limits<double>::infinity();
    
    // Get the vector for the data point
    NumericVector point = data(i, _);
    
    // Check each nueron for bmu
    for (int j = 0; j < neurons.nrow(); j++) {
      // Get the vector for the neuron
      NumericVector neuron = neurons(j, _);
      
      // Calculate the Euclidean distance between the data and neuron vectors
      double dist = calc_dist(neuron, point);
      
      // Check for better distance
      if (dist < bmu_dist)
        bmu_dist = dist;
    } // end for (j)
    
    total_dist = total_dist + bmu_dist;	
  } // end for (i)
  
  // Calculate the error
  double err = total_dist / data.nrow();
  
  // Return list
  List out = List::create(Named("val") = err);
  
  return out;
} // end GetQuantizationError

// Reference:
// G. Polzlbauer, Survey and comparison of quality measures for self-organizing maps,
// in Proc. 5th Workshop Data Analysis, pg 67–82, 2004.
// [[Rcpp::export(name = "get.top.err")]]
List GetTopographicError(List map) {
  // Create matrices for the data and neurons
  NumericMatrix data = df_to_mat(DataFrame(map["data"]));
  NumericMatrix neurons = df_to_mat(DataFrame(map["neurons"]));
    
  // Create variables for other map properties we need
  int x = map["xdim"];
  int N = data.nrow();
    
  // Initialize error count
  int errors = 0;

  // Check each data point
  for (int i = 0; i < N; i++) {
    // Intialize variables
    int bmu_index = 0;
    int sbmu_index = 0;
    double bmu_dist = std::numeric_limits<double>::infinity();
    double sbmu_dist = std::numeric_limits<double>::infinity();
    
    // Get the vector for the data point
    NumericVector point = data(i, _);
    
    // Check each nueron for bmu
    for (int j = 0; j < neurons.nrow(); j++) {
      // Get the vector for the neuron
      NumericVector neuron = neurons(j, _);
      
      // Calculate the euclidean distance between the data and neuron vectors
      double dist = calc_dist(neuron, point);
      
      // Check for better best distance (if better than bmu, its better than sbmu too)
      if (dist < bmu_dist) {
        // Update distances
        sbmu_dist = bmu_dist;
        bmu_dist = dist;
        
        // Update indices
        sbmu_index = bmu_index;
        bmu_index = j;
      } // If not, check second best distance
      else if (dist < sbmu_dist) {
        // Update sbmu distance
        sbmu_dist = dist;
          
        // Update index
        sbmu_index = j;
      } // end if
    } // end for (j)
 
    // Considering the neighborhood:
    // n-xdim-1 n-xdim n-xdim+1
    //   n-1       n      n+1
    // n+xdim+1 n+xdim n+xdim+1
    
    // Find index difference
    int dif = abs(bmu_index - sbmu_index);
    
    // Check for error
    if (!(dif == 1 || dif == x - 1 || dif == x || dif == x + 1))
        errors++;	
  } // end for (i)
  
  // Calculate the error
  double err = (double)errors / N;
  
  // Return list
  List out = List::create(Named("val") = err);
  
  return out;
} // end GetTopographicError

// Reference:
// T.  Villmann, R.  Der, M.  Herrmann, and T.  Martinetz, Topology preservation in self-organizing 
// feature maps: exact definition and measurement, IEEE Trans. Neural Netw., vol. 8 no. 2, 
// pg 256 - 266, 1997.
// [[Rcpp::export(name = "get.top.func")]]
List GetTopographicFunction(List map) {
  // Create matrices for the data and neurons
  NumericMatrix data = df_to_mat(DataFrame(map["data"]));
  NumericMatrix neurons = df_to_mat(DataFrame(map["neurons"]));
  
  // Create variables for other map properties we need
  int x = map["xdim"];
  int y = map["ydim"];
  int N = data.nrow();
  
  // Initialize the connectivity and Delaunay Triangulation matrices
  NumericMatrix C(x * y, x * y);
  NumericMatrix Dm(x * y, x * y);
  
  // Build the connectivity matrix
  for (int i = 0; i < N; i++) {
    // Intialize variables
    int bmu_index = 0;
    int sbmu_index = 0;
    double bmu_dist = std::numeric_limits<double>::infinity();
    double sbmu_dist = std::numeric_limits<double>::infinity();
    
    // Get the vector for the data point
    NumericVector point = data(i, _);
   
    // Check each nueron for bmu
    for (int j = 0; j < neurons.nrow(); j++) {
      // Get the vector for the neuron
      NumericVector neuron = neurons(j, _);
      
      // Calculate the euclidean distance between the data and neuron vectors
      double dist = calc_dist(neuron, point);
      
      // Check for better best distance (if better than bmu, its better than sbmu too)
      if (dist < bmu_dist) {
        // Update distances
        sbmu_dist = bmu_dist;
        bmu_dist = dist;
        
        // Update indices
        sbmu_index = bmu_index;
        bmu_index = j;
      } // If not, check second best distance
      else if (dist < sbmu_dist) {
        // Update sbmu distance
        sbmu_dist = dist;
        
        // Update index
        sbmu_index = j;
      } // end if
      
      //Rcout << "j " << j << std::endl;
    } // end for (j)
    
    C(bmu_index, sbmu_index) = 1;
    C(sbmu_index, bmu_index) = 1;
  } // end for (i)
  
  // Build the Delaunay Triangulation matrix (shortest paths)
  
  // TODO: Need to implement
  
  // Note: Limitation of top func is number of input vectors necessary to compute Dm (see pg 260)
  
  // Initialize function results
  
  NumericVector ks = NumericVector(2 * x * y - 1);
  NumericVector phi = NumericVector(2 * x * y - 1);
  
  //Rcout << "got here" << std::endl;
  
  // Calculates all function values
  for (int i = 0; i < 2 * x * y - 1; i++) {
    int k = i - x * y + 1;
    double p = 0.0;
    
    // Calculate phi(k)
    for (int j = 0; j < x * y; j++) {
      for (int l = 0; l < x * y; l++) {
        int f = 0;
        NumericVector i_idx(2);
        NumericVector j_idx(2);
        i_idx(0) = j % x;
        i_idx(1) = j / x;
        j_idx(0) = l % x;
        j_idx(1) = l / x;
      
        // Calculate f(k)
        if (k > 0) {
          double dist = calc_dist_max(i_idx, j_idx);
          double dist_Dm = 1.0;
          if (dist > k && dist_Dm == 1)
            f++;
        } else if (k < 0) {
          double dist = calc_dist(i_idx, j_idx);
          double dist_Dm = k + 1;
          if (dist == 1 && dist_Dm > k)
            f++;
        } // end if
        
        p = p + f;
      } // end for (l)
    } // end for (j)

    ks(i) = k;
    phi(i) = p / (x * y);
  } // end for (i)
  
  // Set phi(0)
  phi(x * y) = phi(x * y + 1) + phi(x * y - 1); 
  
  // Return list
  List out = List::create(Named("C") = C,
                          Named("k") = ks,
                          Named("phi") = phi);
  
  return out;
} // end GetTopographicFunction

// Reference:
// J. Venna and S. Kaski, "Neighborhood preservation in nonlinear projection methods: 
// An experimental study", Lecture Notes in Comput. Sci., vol. 2130, pg 485-491, 2001.
// [[Rcpp::export(name = "get.hood.pres")]]
List GetNeighborhoodPreservation(List map, int k) {
  // Initialize
  double M_1 = 0.0;
  double M_2 = 0.0;
  
  // Create matrices for the data and neurons
  NumericMatrix data = df_to_mat(DataFrame(map["data"]));
  NumericMatrix neurons = df_to_mat(DataFrame(map["neurons"]));
  NumericVector projection = df_to_mat(DataFrame(map["visual"]));
  
  // Create list of code book vectors from the projection
  NumericMatrix pdata(projection.size(), neurons.cols());
  for (int i = 0; i < projection.size(); i++) {
    int p = projection[i] - 1;
    pdata(i, _) = neurons(p, _);
  } // end for (i)
  
  // Get the distances for the data and projection
  NumericMatrix dist_mat = calc_dist_mat(data);
  NumericMatrix pdist_mat = calc_dist_mat(pdata);
  
  // Check each data point
  for (int i = 0; i < data.nrow(); i++) {
    // Get the distances for x_i
    NumericVector dist = dist_mat.row(i);
    NumericVector pdist = pdist_mat.row(i);
    
    // Initialize the unsorted index vectors
    std::vector<int> idx(dist.size());
    std::vector<int> pidx(pdist.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::iota(pidx.begin(), pidx.end(), 0);
    
    // Sort the indices by the distance
    std::sort(idx.begin(), idx.end(),
                [dist](double i1, double i2) {return dist[i1] < dist[i2];});
    std::sort(pidx.begin(), pidx.end(),
                [pdist](double i1, double i2) {return pdist[i1] < pdist[i2];});
    
    // Get all x_j in (and not in)  C_k(x_i) 
    std::vector<int> Ck(idx.begin(), idx.begin() + k);
    std::vector<int> not_Ck(idx.begin() + k, idx.end());
    
    // Get all x_j in (and not in) C^_k(x_i)
    std::vector<int> hat_Ck(pidx.begin(), pidx.begin() + k);
    std::vector<int> not_hat_Ck(pidx.begin() + k, pidx.end());
    
    // Get U_k(x_i), e.g. the intersection of x_j not in C_k(x_i) and x_j in C^_k(x_i)
    std::vector<int> Uk(k);
    std::sort(not_Ck.begin(), not_Ck.end());
    std::sort(hat_Ck.begin(), hat_Ck.end());
    std::vector<int>::iterator it = std::set_intersection(not_Ck.begin(), not_Ck.end(), 
                                                          hat_Ck.begin(), hat_Ck.end(), Uk.begin());
    Uk.resize(it - Uk.begin());
    
    // Get V_k(x_i), e.g. the intersection of x_j in C_k(x_i) and x_j not in C^_k(x_i)
    std::vector<int> Vk(k);
    std::sort(Ck.begin(), Ck.end());
    std::sort(not_hat_Ck.begin(), not_hat_Ck.end());
    it = std::set_intersection(Ck.begin(), Ck.end(), not_hat_Ck.begin(), not_hat_Ck.end(), Vk.begin());
    Vk.resize(it - Vk.begin());
    
    // Calculate the inner sums
    for (int j = 0; j < std::max(Uk.size(), Vk.size()); j++) {
      if (j < Uk.size()) {
          int x_j = Uk[j];
          if (x_j != i) {
            int r = find(idx.begin(), idx.end(), x_j) - idx.begin() + 1;
            M_1 = M_1 + r - k;
          } // end if (x_j)
      } // end if (j)
      if (j < Vk.size()) {
        int x_j = Vk[j];
        if (x_j != i) {
          int r_hat = find(pidx.begin(), pidx.end(), x_j) - pidx.begin() + 1;
          M_2 = M_2 + r_hat - k;
        } // end if (x_j)
      } // end if (j)
    } // end for (j)
  } // end for (i)
  
  // Convert the sum
  int N = data.nrow();
  M_1 = 1 - (2 * M_1 / (N * k * (2 * N - 3 * k - 1)));
  M_2 = 1 - (2 * M_2 / (N * k * (2 * N - 3 * k - 1)));
  
  // Return list
  List out = List::create(Named("k") = k,
                          Named("trustworthiness") = M_1,
                          Named("neighborhood.preservation") = M_2);
  
  return out;
} // end GetNeighborhoodPreservation

// Reference:
// J. Hirschberg and A. Rosenberg, V-Measure: A conditional entropy-based external cluster evaluation,
// Columbia University Academic Commons, 2007, http://hdl.handle.net/10022/AC:P:21139
// [[Rcpp::export(name = "get.v.measure")]]
List GetVMeasure(IntegerVector labels, IntegerVector clusters, double beta = 1.0) {
  // Check the sizes
  if (labels.size() != clusters.size())
    stop("get.v.measure: label and cluster vectors sizes don't match.");
  
  // Get the levels
  IntegerVector label_lev = sort_unique(labels);
  IntegerVector clust_lev = sort_unique(clusters);
  
  // Initialize sizes
  int N = labels.size();
  int n = label_lev.size();
  int m = clust_lev.size(); 
  
  // Generate the contingency table
  NumericMatrix A(n, m);
  for (int i = 0; i < N; i++)
    A(labels[i] - 1, clusters[i] - 1) = A(labels[i] - 1, clusters[i] - 1) + 1;
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
        H_CK = H_CK + A(c, k) * (log(A(c, k)) - log(sum(A(_, k))));
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
        H_KC = H_KC + A(c, k) * (log(A(c, k)) - log(sum(A(c, _))));
  H_KC = - H_KC;
  
  // Calculate H(K)
  for (int k = 0; k < m; k++)
    H_K = H_K + sum(A(_, k)) * log(sum(A(_, k)));
  H_K = - H_K;
  if (std::isnan(H_K))
    H_K = 0;
  
  // Calculate homogeneity
  if (H_C == 0)
    homo = 1;
  else
    homo = 1 - H_CK / H_C;
  
  // Calculate completeness
  if (H_K == 0)
    comp = 1;
  else
    comp = 1 - H_KC / H_K;
  
  // Calculate the weighted harmonic mean of homogeneity and completeness
  double v = ((1 + beta) * homo * comp) / ((beta * homo) + comp);
  
  // Return list
  List out = List::create(Named("beta") = beta,
                          Named("homogeneity") = homo,
                          Named("completeness") = comp,
                          Named("v.measure") = v);
  
  return out;
} // end GetVMeasure
