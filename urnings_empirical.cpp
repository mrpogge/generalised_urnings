#include <Rcpp.h>
#include <algorithm>  // std::upper_bound
#include <vector>
#include <cmath>      // std::floor
using namespace Rcpp;

// [[Rcpp::export]]
List urnings_empirical(IntegerVector X,// item scores
                       IntegerVector Pers,// person ID [should start with 0]
                       IntegerVector Item,// item ID [should start with 0]
                       int n_students,//P
                       int n_items,//I
                       int n_pars, //J
                       int n_dims, //D
                       int n_resp, //total number of response in the dataset
                       int n_categories, // M
                       int n_scores, //K
                       IntegerVector start_u, // starting values for the student urnings (flattened PxD matrix)
                       IntegerVector start_v,// starting values for the item urnings (flattened Ix(J-D) matrix)
                       IntegerVector S, // flattened scoring matrix connecting categories and parameters
                       IntegerVector Smax, // vector with column-wise maxima of matrix s
                       IntegerVector urn_size,// J-length vector 
                       IntegerVector blocks, // column IDs of the columns of the combined b-matrix for where the details belonging to each of the scores in tilde_t begin
                       IntegerVector Upd,// matrices b put next to each other and flattened into a vector
                       IntegerVector Score, // vector with all possible non-zero updates (tilde_t)
                       IntegerVector paired) {// for each item parameter specifies whether paired update is done for it (where these parameter is an element of w)
  
  GetRNGstate(); // Access R's RNG
  int n_itempars=n_pars-n_dims;
  
  // creating vectors which will be returned to R
  // urnings before the update 
  IntegerVector RR(n_resp*n_pars);
  // urnings after the update
  IntegerVector RR_after(n_resp*n_pars);
  // person-specific MSE of the item scores
  NumericVector MSE_p_score(n_students,0.0);
  // person-specific MSE of the binary scores
  NumericVector MSE_p_binary(n_students,0.0);
                                         
  // Expected values for binary score for each response
  NumericVector Expected(n_resp,0.0);
  
  // vector to count the number of response each student had, needed to compute their MSE
  std::vector<int> student_n_resp(n_students,0);
  
  // vectors for internal C computations
  
  // student urnings; filled up with the starting values
  std::vector<int> u(n_students * n_dims);
  for (int k = 0; k < n_dims; ++k) {
    int dim_offset = k * n_students;
    for (int i = 0; i < n_students; ++i) {
      u[dim_offset + i] = start_u[k];
    }
  }
  
  // item urnings; filled up with the starting values
  std::vector<int> v(n_items * n_itempars);
  for (int k = 0; k < n_itempars; ++k) {
    int par_offset = k * n_items;
    for (int i = 0; i < n_items; ++i) {
      v[par_offset + i] = start_v[k];
    }
  }
  
  // duplicate R input with C-vectors (saves computation time)
  std::vector<int> urnSize(urn_size.begin(), urn_size.end());
  std::vector<int> Blocks(blocks.begin(), blocks.end());
  std::vector<int> upd(Upd.begin(), Upd.end());
  std::vector<int> score(Score.begin(), Score.end());
  std::vector<int> Paired(paired.begin(), paired.end());
  
  //vectors to save CDFs for sampling items and responses.
  std::vector<double> cmsm(n_items + 1, 0.0);
  std::vector<double> cmsm2(n_categories+1, 0.0);
  cmsm2[n_categories]=1.0;
  
  std::vector<int> idx_i(n_dims,0);
  std::vector<int> idx_j(n_itempars,0);
  std::vector<int> C(n_items, 0);
  
  // vector to store current ratings
  std::vector<int> R(n_pars, 0);
  
  // creating the queues for the items 
  std::vector<int> queue(n_items * n_itempars, 0);
  
  std::vector<int> SW(n_items * n_scores);
  std::vector<int> SWsum(n_scores);
  
  std::vector<int> sel(n_pars+1,0);
  sel[0]=1;
  for(int i=0;i<n_pars;++i){
    sel[i+1]=sel[i]*(urnSize[i]+1);
  }
  
  std::vector<int> sel2(n_pars+1,0);
  sel2[0]=1;
  for(int i=0;i<n_pars;++i){
    sel2[i+1]=sel2[i]*(urnSize[i]+Smax[i]+1);
  }
  
  double p=0.0;
  int idx_pers=0, idx_item=0, idx=0;
  int x=0, y=0, j=0;
  int accept=0;
  
  double MSE_binary=0.0;
  double MSE_score=0.0;  
  
  // loop over all responses
  for(int g=0;g<n_resp;++g){
    int i=Pers[g];
    int j=Item[g];
    // count the number of responses of the student
    student_n_resp[i]+=1;
    
    // save the indices for accessing the current ratings of person i in replication r in dimentions from 0 to n_dims-1
    for(int k=0;k<n_dims;++k){
      idx_i[k]=i+k*n_students;
    }
    
    // collect current ratings of the person and of the item in vector R. The old values are left in u and v
    for (int k=0; k<n_dims; ++k) {
      R[k]=u[idx_i[k]];
    }
    for (int k=n_dims; k<n_pars; ++k) {
      int w=k-n_dims;
      idx_j[w]=j+w*n_items;
      R[k]=v[idx_j[w]];
    }
    
    // compute probabilities for expected scores using the urnings before the update
    NumericVector P_before(n_categories,0.0);
    double  CP=0.0;
    for(int m=0;m<n_categories;++m){
      double  PP=1.00;
      for(int k=0;k<n_pars;++k){
        int r=R[k];
        int smax=Smax[k];
        int n=urn_size[k];
        int s=S[m+n_categories*k];
        if(s>0){
          for(int h=0;h<s;++h){
            PP*=(r-h);
          }
        }
        if(s<smax){
          for(int h=0;h<(smax-s);++h){
            PP*=n-r-h;
          }
        }
      }
      P_before[m]=PP;
      CP+=PP;
    }
    for(int m=0;m<n_categories;++m){// normalise the probabilities
      if((P_before[m]==0.0)&(CP==0.0)){P_before[m]=1.0/n_categories;}else{
        P_before[m]=P_before[m]/CP;
      }
    }
    
    // compute the expected score
    double expected_score=0.0;
    for(int m=1;m<n_categories;++m){
      expected_score+=P_before[m]*m;
    }
    // compute expected binary score
    double expected_binary_score=P_before[n_categories-1];
    if(n_categories>2){expected_binary_score+=P_before[n_categories-2];}
    
    // observed score
    x=X[g];
    
    double error=expected_score-x;
    double error2=error;
    if(n_categories>2){
      error2=expected_binary_score-1*(x>(n_categories-3));
    }
    Expected[g]=expected_binary_score;
    MSE_p_score[i]+=error*error;
    MSE_score+=error*error;
    MSE_p_binary[i]+=error2*error2;
    MSE_binary+=error2*error2;
    
    //save the values before the update	
    for(int k=0;k<n_pars;++k){
      RR[k+g*n_pars]=R[k];
    }
    
    // Add x-scores to ratings
    for (int k=0; k<n_pars; ++k) {
      R[k]+=S[x+n_categories*k];
    }
    
    // compute probabilities to simulate
    NumericVector P(n_categories,0.0);
    CP=0.0;
    for(int m=0;m<n_categories;++m){
      double  PP=1.00;
      for(int k=0;k<n_pars;++k){
        int r=R[k];
        int smax=Smax[k];
        int n=urn_size[k]+smax;
        int s=S[m+n_categories*k];
        if(s>0){
          for(int h=0;h<s;++h){
            PP*=(r-h);
          }
        }
        if(s<smax){
          for(int h=0;h<(smax-s);++h){
            PP*=n-r-h;
          }
        }
      }
      P[m]=PP;
      CP+=PP;
    }
    
    // normalize probabilities
    for(int m=0; m<(n_categories-1); ++m) {
      if((P[m] == 0.0) & (CP == 0.0)) {
        P[m] = 1.0/n_categories;
      } else {
        P[m] = P[m]/CP;
      }
    }
    
    // construct cumulative probabilities
    cmsm2[0] = 0.0;
    
    for(int m=0; m<(n_categories-1); ++m) {
      cmsm2[m+1] = cmsm2[m] + P[m];
    }
    
    // Sample sim response 
    p=R::runif(0.0,1.0);
    y=0;
    for (int s=0; s<n_categories; ++s) {
      if (p <= cmsm2[s+1]) {
        y=s;
        break;
      }
    }
    
    // Subtract y-scores from the ratings
    for (int k=0; k<n_pars; ++k) {
      R[k]-=S[y+n_categories*k];
    }
    
    for(int k=0;k<n_dims;k++){// update the person ratings
      u[idx_i[k]]=R[k];
    }
    // update the item rating taking paired update into account
    for (int w = 0; w < n_itempars; ++w) {// for the item parameters it depends on whether they are in Paired update
      if(Paired[w]==0){
        v[idx_j[w]]=R[w+n_dims];// if the item parameter is not part of the update scheme, then update it directly
      }else{ // if they are in Paired update, it depends on whether there are items waiting for an update in a different direction
        int Old=v[idx_j[w]];
        int New=R[w+n_dims];
        if (Old!=New){ // this is needed because a specific item parameter might have a zero update even if x!=y
          int change = New - Old;
          int s = 0;
          for (int i = 0; i < n_scores; ++i) {// which element of score is change?
            if (score[i] == change) {
              s = i;
              break;
            }
          }
          
          // Reset SW and SWsum (sampling weights for selecting an item from the queue) 
          std::fill(SW.begin(), SW.end(), 0); 
          std::fill(SWsum.begin(), SWsum.end(), 0);
          
          // Precompute constants
          int base_idx = w * n_items;
          int n = urnSize[n_dims + w];
          
          // Fill SW and SWsum
          for (int t = 0; t < n_scores; ++t) {
            int score_t = score[t];
            int offset = t * n_items;
            
            for (int ii = 0; ii < n_items; ++ii) {
              idx = ii + base_idx;
              int q = queue[idx];
              int vv = v[idx];
              
              // Compute validity as an int (0 or 1)
              int valid = (q == score_t) &&
                ((vv + score_t) <= n) &&
                ((vv + score_t) > -1) &&
                (ii != j);
              
              SW[ii+offset] = valid;
              SWsum[t] += valid;
            }
          }
          
          
          int success = 0;
          for (int d = Blocks[s]; d < Blocks[s + 1]; ++d) {// go through all possible configurations for the Paired update when the update was s
            if (success == 1) break;
            
            bool enough = true;
            for (int t = 0; t < n_scores; ++t) {
              if (SWsum[t] < upd[t+d*n_scores]) {
                enough = false;
                break;
              }
            }
            
            if (enough) {
              success = 1;
              for (int t = 0; t < n_scores; ++t) {
                int count = upd[t+d*n_scores];
                if (count > 0) {// if update t is needed
                  // Build cumsum for sampling items from the queue
                  cmsm[0] = 0;
                  for (int jj = 0; jj < n_items; ++jj) {
                    cmsm[jj + 1] = cmsm[jj] + static_cast<double>(SW[jj + t * n_items]);
                  }
                  
                  for(int h=0;h<count;++h){
                    p = R::runif(0, cmsm[n_items]);// sample one item from the queue
                    int jj = std::upper_bound(cmsm.begin(), cmsm.end(), p) - cmsm.begin() - 1;
                    
                    idx = jj + w * n_items;
                    v[idx] += score[t];// update the sampled item
                    queue[idx] = 0;// remove from the queue
                    
                    // update samling weights and the cumulative probabilities 
                    SWsum[t] -= 1;
                    SW[jj + t * n_items] = 0;
                    for (int k = jj + 1; k <= n_items; ++k) {
                      cmsm[k] -= 1.0;
                    }
                  }	  
                }
              }
              v[idx_j[w]] += score[s];// the item parameter is updted if we should a match in the queue
            }
          }
          if (success == 0) {
            queue[idx_j[w]] = score[s];// if we did not fill a match, no update happens and the item is put in the queue itself
          }
        }  
      }
    }
    
    for(int k=0;k<n_pars;++k){
      RR_after[k+g*n_pars]=R[k];//save the values before the update	
    }
    
  }
  
  for(int i=0;i<n_students;i++){
    MSE_p_score[i]/=student_n_resp[i];
    MSE_p_binary[i]/=student_n_resp[i];
  }
  MSE_score/=n_resp;
  MSE_binary/=n_resp;
  
  PutRNGstate(); // Save RNG state
  
  return Rcpp::List::create(
    Rcpp::Named("Ratings_before") = RR,
    Rcpp::Named("Ratings_after") = RR_after,
    Rcpp::Named("MSE_p_score")=MSE_p_score,
    Rcpp::Named("MSE_p_binary")=MSE_p_binary,
    Rcpp::Named("MSE_score")=MSE_score,
    Rcpp::Named("MSE_binary")=MSE_binary,
    Rcpp::Named("Expected")=Expected
  );
}    