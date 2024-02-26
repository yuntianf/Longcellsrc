#include "umi_dist.h"
#include "reads_filter.h"



int iso2_dis(std::string a,std::string b,
             const std::string split,
             const std::string sep){
  if(a == b){
    return(0);
  }
  std::vector<int> sites_a = isoform2sites(a, split,sep);
  std::vector<int> sites_b = isoform2sites(b, split,sep);

  int a_size = sites_a.size();
  int b_size = sites_b.size();

  int i =0,j = 0, intersect = 0;
  while(i < a_size/2){
    if(sites_a[2*i] > sites_b[2*j+1]){
      if(j == b_size/2-1){
        break;
      }
      j++;
    }
    else if(sites_a[2*i+1] < sites_b[2*j]){
      i++;
    }
    else{
      intersect += bin2_intersect(sites_a[2*i],sites_a[2*i+1],
                                  sites_b[2*j],sites_b[2*j+1]);
      i++;
    }
  }

  int diff = iso_len(sites_a) + iso_len(sites_b) - 2*intersect;
  return(diff);
}
// [[Rcpp::export]]
DataFrame isos_dis(const std::vector<std::string> isoforms,const int thresh,
             const std::string split,const std::string sep){
  int n = isoforms.size();

  std::vector<int> node1;
  std::vector<int> node2;
  std::vector<int> dis;
  for(int i =0;i < n;i++){
    for(int j = i;j < n;j++){
      int temp = iso2_dis(isoforms[i],isoforms[j],split,sep);
      if(temp <= thresh){
        node1.push_back(i+1);
        node2.push_back(j+1);
        dis.push_back(temp);
      }
    }
  }

  DataFrame out = DataFrame::create(Named("node1") = node1 , Named("node2") = node2, Named("dis") = dis);
  return(out);
}

// [[Rcpp::export]]
NumericVector size_filter_cpp(NumericVector size,double ratio){
  int n = size.size();
  double left = 0, right = sum(size);
  int id = 0;
  double diff = 100;
  for(int i = 0;i < n;i++){
    left += size[i];
    right -= size[i];
    double thresh = left-right*ratio/(1-ratio);
    //Rcout << left<< "-" << right*ratio/(1-ratio) << endl;
    if(diff > abs(thresh)){
      diff = abs(thresh);
      id = i;
    }
    if(thresh >= 0){
      break;
    }
  }
  NumericVector weight(n,0.0);
  weight[size > size[id]] = 1;
  weight[size == size[id]] = (sum(size <= size[id])-double(id)-1.0)/sum(size == size[id]);
  return(weight);
}
