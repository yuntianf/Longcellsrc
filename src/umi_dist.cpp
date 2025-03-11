#include "umi_dist.h"
#include "edit.h"

int UMI_NS[20000][20000] = {};
int UMI_CORRES[400000] = {};
static int none = -100;

// ##### isoform cluster #####
std::vector<std::string> flatten(std::vector<std::vector<std::string> > const &vec){
  std::vector<std::string> flattened;
  for (auto const &v: vec) {
    flattened.insert(flattened.end(), v.begin(), v.end());
  }
  return flattened;
}


std::vector<std::string> str_split(std::string s, std::string split){
  std::vector<std::string> out;

  size_t pos = 0;
  std::string token;
  while ((pos = s.find(split)) != std::string::npos) {
    token = s.substr(0, pos);
    out.push_back(token);
    s.erase(0, pos + split.length());
  }
  if(s.size() > 0){
    out.push_back(s);
  }

  return(out);
}

// [[Rcpp::export]]
std::vector<int> isoform2sites(std::string iso,
                                      const std::string split,
                                      const std::string sep){
  std::vector<std::string> exons = str_split(iso,split);

  std::vector<std::vector<std::string> > bins;
  int len = exons.size();
  for(int i = 0;i < len;i++){
    bins.push_back(str_split(exons[i],sep));
  }
  std::vector<std::string> sites_chr = flatten(bins);

  std::vector<int> sites_int;
  transform(sites_chr.begin(), sites_chr.end(), back_inserter(sites_int),
            [](const std::string& str) { return stoi(str); });

  return(sites_int);
}


int iso_len(std::vector<int> sites){
  int sites_size = sites.size();

  if(sites_size % 2 > 0){
    cout << "the length of sites should be an even" << endl;
    exit(1);
  }
  int len = 0;
  for(int i = 0;i < sites_size/2;i++){
    int sites_s = sites[i*2];
    int sites_e = sites[i*2 + 1];

    len= len+sites_e-sites_s+1;
  }
  return(len);
}


int bin2_intersect(int a_start,int a_end,
                          int b_start,int b_end){
  int len = 0;
  if(a_start <= b_end && a_end >= b_start){
    len = min({a_end,b_end})-max({a_start,b_start}) + 1;
  }
  return(len);
}

// [[Rcpp::export]]
std::vector<int> sites_chop(std::vector<int> sites, int start,int end){
  int sites_size = sites.size();

  if(sites_size % 2 > 0){
    cout << "the length of sites should be an even" << endl;
    exit(1);
  }

  std::vector<int> chop;
  for(int i = 0;i < sites_size/2;i++){
    int sites_s = sites[i*2];
    int sites_e = sites[i*2 + 1];
    if(sites_e <= start || sites_s >= end){
      continue;
    }
    else{
      int temp = max({sites_s,start});
      chop.push_back(max({sites_s,start}));
      chop.push_back(min({sites_e,end}));
    }
  }
  return(chop);
}

// [[Rcpp::export]]
int iso2_mid_dist(std::string a,std::string b,
                    const std::string split,
                    const std::string sep){
  if(a == b){
    return(0);
  }
  std::vector<int> sites_a = isoform2sites(a, split,sep);
  std::vector<int> sites_b = isoform2sites(b, split,sep);

  int a_size = sites_a.size();
  int b_size = sites_b.size();

  int start = max({sites_a[0],sites_b[0]});
  int end = min({sites_a[a_size - 1],sites_b[b_size - 1]});

  if(start >= end){
    return(-1);
  }

  std::vector<int> chop_a = sites_chop(sites_a, start, end);
  std::vector<int> chop_b = sites_chop(sites_b, start, end);

  int chop_a_size = chop_a.size();
  int chop_b_size = chop_b.size();

  int intersect = 0;
  int next_start = 0;
  for(int i =0;i < chop_a_size/2;i++){
    int j = next_start;
    bool flag = false;

    while(j < chop_b_size/2){
      int a_start = chop_a[2*i];
      int a_end = chop_a[2*i + 1];
      int b_start = chop_b[2*j];
      int b_end = chop_b[2*j + 1];

      int inter = bin2_intersect(a_start,a_end,b_start,b_end);

      if(inter > 0){
        if(!flag){
          next_start = j;
          flag = true;
        }
        intersect += inter;
      }
      else{
        if(flag){
          break;
        }
      }
      j++;
    }
  }

  int dis = iso_len(chop_a) + iso_len(chop_b) - 2*intersect;
  return(dis);
}

NumericMatrix iso_mid_dist(std::vector<std::string> iso_set,
                           int thresh,
                           std::string split,std::string sep){
  int iso_size = iso_set.size();

  NumericMatrix iso_dis ((iso_size+1)*iso_size/2,2);
  int id = 0;

  for(int i = 0;i < iso_size;i++){
    for(int j = i;j < iso_size;j++){
      if(i == j){
        iso_dis(id,0) = i;
        iso_dis(id,1) = j;
        id = id + 1;
      }
      else{
        int dis = iso2_mid_dist(iso_set[i],iso_set[j],
                           split,sep);
        //cout << i << "-" << j<<":" << dis << endl;
        if(dis >=0 && dis <= thresh){
          iso_dis(id,0) = i;
          iso_dis(id,1) = j;
          id = id+ 1;
        }
      }
    }
  }
  return(iso_dis(Range(0,id-1),_));
}

// [[Rcpp::export]]
NumericVector iso2_mid_diff(std::string a,std::string b,
                            const int end_bias,
                 const std::string split,
                 const std::string sep){
  NumericVector out(2);
  if(a == b){
    out = {0,1};
    return(out);
  }
  std::vector<int> sites_a = isoform2sites(a, split,sep);
  std::vector<int> sites_b = isoform2sites(b, split,sep);

  int a_size = sites_a.size();
  int start = sites_a[0];
  int end = sites_a[a_size - 1];

  std::vector<int> chop_b = sites_chop(sites_b, start, end);

  int b_size = chop_b.size();

  if(b_size == 0){
    out = {-1,-1};
    return(out);
  }

  if(start <= chop_b[0] && start >= chop_b[0]-end_bias &&
    chop_b[0] <= sites_a[1]){
    start = chop_b[0];
  }
  if(end >= chop_b[b_size - 1] && end <= chop_b[b_size - 1]+end_bias &&
     chop_b[b_size - 1] >= sites_a[a_size - 2]){
    end = chop_b[b_size - 1];
  }

  if(start >= end){
    out = {-1,-1};
    return(out);
  }


  //Rcout << start << " " << end << endl;

  std::vector<int> chop_a = sites_chop(sites_a, start, end);

  int chop_a_size = chop_a.size();
  int chop_b_size = chop_b.size();

  /*
  for (int i = 0; i < chop_a_size; ++i) {
    Rcpp::Rcout << chop_a[i] << " ";
  }
  Rcout<< endl;
  for (int i = 0; i < chop_b_size; ++i) {
    Rcpp::Rcout << chop_b[i] << " ";
  }
  Rcout<< endl;
   */

  int intersect = 0;
  int next_start = 0;
  for(int i =0;i < chop_a_size/2;i++){
    int j = next_start;
    bool flag = false;

    while(j < chop_b_size/2){
      int a_start = chop_a[2*i];
      int a_end = chop_a[2*i + 1];
      int b_start = chop_b[2*j];
      int b_end = chop_b[2*j + 1];

      int inter = bin2_intersect(a_start,a_end,b_start,b_end);

      if(inter > 0){
        if(!flag){
          next_start = j;
          flag = true;
        }
        intersect += inter;
      }
      else{
        if(flag){
          break;
        }
      }
      j++;
    }
  }

  // cout << iso_len(chop_a) << " " << iso_len(chop_b) << " " << intersect << endl;

  double dis = iso_len(chop_a) + iso_len(chop_b) - 2*intersect;
  double ratio = double(intersect)/double(iso_len(sites_b));

  out = {dis,ratio};
  return(out);
}

// [[Rcpp::export]]
NumericMatrix isoset_mid_diff(std::vector<std::string> iso_set1,
                              std::vector<std::string> iso_set2,
                              const int thresh,
                              const double overlap_thresh,
                              const int end_bias,
                              std::string split,
                              std::string sep){
  int iso1_size = iso_set1.size();
  int iso2_size = iso_set2.size();

  NumericMatrix iso_dis (iso1_size*iso2_size,4);
  int id = 0;

  for(int i = 0;i < iso1_size;i++){
    for(int j = 0;j < iso2_size;j++){
      NumericVector out = iso2_mid_diff(iso_set1[i],iso_set2[j],end_bias,
                               split,sep);
      int dis = out[0];
      double ratio = out[1];
        //cout << i << "-" << j<<":" << dis << endl;
      if(dis >=0 && dis <= thresh && ratio >= overlap_thresh){
        iso_dis(id,0) = double(i);
        iso_dis(id,1) = double(j);
        iso_dis(id,2) = double(dis);
        iso_dis(id,3) = double(ratio);
        id = id+ 1;
      }
    }
  }
  if(id == 0){
    iso_dis(id,0) = -1;
    id = 1;
  }
  return(iso_dis(Range(0,id-1),_));
}


// ##### UMI distance #####
// ### needleman ###
int needle(std::string A,std::string B,
                        int match_score,
                        int mismatch_score,
                        int gap_score){
  int n = A.size();
  int m = B.size();

  NumericMatrix dp(n+1,m+1);
  for (int i=0;i<=n;i++) dp(i,0) = dp(0,i) = i * gap_score;
  for (int i=1;i<=n;i++)
  {
    for (int j=1;j<=m;j++)
    {
      int S = (A[i-1] == B[j-1]) ? match_score : mismatch_score;
      dp(i,j) = max({dp(i-1,j-1) + S, max({dp(i-1,j) + gap_score, dp(i,j-1) + gap_score})});
    }
  }

  return dp(n,m);
}

// ### edit distance ###
int minEditDis(std::string seq1, std::string seq2, int k) {
  std::vector<std::string> seq1_kmer = kmer(seq1, k,1);
  std::vector<std::string> seq2_kmer = kmer(seq2, k,1);

  int kmer_size = seq1_kmer.size();

  int edit = seq1.size();
  for(int i = 0;i < kmer_size;i++){
    for(int j = 0;j < kmer_size;j++){
      int dis = editDist(seq1_kmer[i],seq2_kmer[j]);
      if(edit > dis){
        edit = dis;
      }
      if(dis == 0){
        break;
      }
    }
  }
  return(edit);
}


// ### index ###
int pair2id(int x, int y){
  int m = max({x,y});
  int n = min({x,y});
  int loc = (m+1)*m/2;
  int id = loc+n;
  return(id);
}

List index(std::vector<std::string> data, std::vector<std::string> uniq){
  List L;
  int uniq_size = uniq.size();
  int data_size = data.size();

  for(int i =0;i < uniq_size;i++){
    std::vector<int> sub_index;
    for(int j = 0;j <data_size;j++){
      if(data[j] == uniq[i]){
        sub_index.push_back(j);
      }
    }
    L.push_back(sub_index);
  }

  return(L);
}

std::vector<std::string> vec_extract(std::vector<std::string> data, std::vector<int> index){
  int len = index.size();

  std::vector<std::string> extract;
  for(int i = 0;i < len;i++){
    extract.push_back(data[index[i]]);
  }

  return(extract);
}


std::vector<std::vector<int> > umi_needle(std::vector<std::string> umi, std::vector<int> count,
                                std::vector<int> umi_id1, std::vector<int> umi_id2,
                                int thresh) {
  int umi_set1_size = umi_id1.size();
  int umi_set2_size = umi_id2.size();

  std::vector<std::vector<int> > umi_ns;

  for (int i = 0; i < umi_set1_size; i++) {
    for (int j = 0; j < umi_set2_size; j++) {
      int ns = none;
      if(UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] != none){
        ns = UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]];
      }
      else {
        ns = needle(umi[umi_id1[i]], umi[umi_id2[j]]);
        UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] = ns;
        UMI_NS[UMI_CORRES[umi_id2[j]]][UMI_CORRES[umi_id1[i]]] = ns;
      }
      if (ns >= thresh) {
        int score = 1;
        if(umi_id1[i] == umi_id2[j]){
          score = count[umi_id1[i]]*(count[umi_id1[i]]-1)/2;
        }
        else{
          score = count[umi_id1[i]]*count[umi_id2[j]];
        }
        std::vector<int> sub_ns = { umi_id1[i],umi_id2[j],ns,score };
        umi_ns.push_back(sub_ns);
      }
    }
  }
  return(umi_ns);
}

std::vector<std::vector<int> > umi_edit(std::vector<std::string> umi,
                                std::vector<int> umi_id1, std::vector<int> umi_id2,
                                int thresh, int k) {
  int umi_set1_size = umi_id1.size();
  int umi_set2_size = umi_id2.size();

  std::vector<std::vector<int> > umi_ns;

  for (int i = 0; i < umi_set1_size; i++) {
    for (int j = 0; j < umi_set2_size; j++) {
      int ns = none;
      if(UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] != none){
        ns = UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]];
      }
      else {
        ns = minEditDis(umi[umi_id1[i]], umi[umi_id2[j]],k);
        UMI_NS[UMI_CORRES[umi_id1[i]]][UMI_CORRES[umi_id2[j]]] = ns;
        UMI_NS[UMI_CORRES[umi_id2[j]]][UMI_CORRES[umi_id1[i]]] = ns;
      }
      if (ns <= thresh) {
        std::vector<int> sub_ns = { umi_id1[i],umi_id2[j],ns };
        umi_ns.push_back(sub_ns);
      }
    }
  }
  return(umi_ns);
}


void initialize(std::vector<std::string> umi){
  std::set<std::string> umi_set(umi.begin(), umi.end());
  std::vector<std::string> umi_uniq;
  umi_uniq.assign(umi_set.begin(), umi_set.end());

  int umi_size = umi.size();
  int umi_uniq_size = umi_uniq.size();
  int umi_len = umi[0].size();

  for(int i = 0;i < umi_uniq_size;i++){
    for(int j = 0;j < umi_uniq_size;j++){
      if(i == j){
        UMI_NS[i][j] = umi_len;
      }
      else{
        UMI_NS[i][j] = none;
      }
    }
  }

  for(int i = 0;i < umi_size;i++){
    for(int j = 0;j < umi_uniq_size;j++){
      if(umi[i] == umi_uniq[j]){
        UMI_CORRES[i] = j;
      }
    }
  }
}


// [[Rcpp::export]]
std::vector<std::vector<int> > umi_graph_table(std::vector<std::string> umi,
                              std::vector<std::string> isoform,
                              std::vector<int> count,
                              int sim_thresh,int iso_thresh,
                              std::string split,std::string sep){
  int umi_len = umi.size();
  int isolen = isoform.size();
  int count_len = count.size();

  if(isolen != umi_len || isolen != count_len || umi_len != count_len){
    cout << "The size of isoforms and umi don't match!" << endl;
  }

  std::set<std::string> iso_set(isoform.begin(), isoform.end());
  std::vector<std::string> iso_uniq;
  iso_uniq.assign(iso_set.begin(), iso_set.end());

  NumericMatrix iso_src = iso_mid_dist(iso_uniq,iso_thresh,
                                   split = split,sep = sep);
  //return(iso_src);
  List iso_index = index(isoform,iso_uniq);
  initialize(umi);

  int iso_pair_count = iso_src.nrow();
  std::vector<std::vector<int> > umi_ns;
  for(int i = 0;i < iso_pair_count;i++){
    std::vector<std::vector<int> > temp_ns = umi_needle(umi,count,
                                       iso_index[iso_src(i,0)],
                                       iso_index[iso_src(i,1)],
                                       sim_thresh);

    if(temp_ns.size() > 0){
      umi_ns.insert(umi_ns.end(),temp_ns.begin(),temp_ns.end());
    }
  }
  return(umi_ns);
}

/*
// [[Rcpp::export]]
std::vector<std::vector<int> > umi_edit_table(std::vector<std::string> umi,
                                     std::vector<std::string> isoform,
                                     int edit_thresh,
                                     int k,
                                     int iso_thresh,
                                     std::string split,std::string sep){
  int umi_len = umi.size();
  int isolen = isoform.size();

  if(isolen != umi_len){
    cout << "The size of isoforms and umi don't match!" << endl;
  }

  set<std::string> iso_set(isoform.begin(), isoform.end());
  std::vector<std::string> iso_uniq;
  iso_uniq.assign(iso_set.begin(), iso_set.end());

  NumericMatrix iso_src = iso_mid_dist(iso_uniq,iso_thresh,
                                   split = split,sep = sep);
  //return(iso_src);
  List iso_index = index(isoform,iso_uniq);
  initialize(umi);

  int iso_pair_count = iso_src.nrow();
  std::vector<std::vector<int> > umi_ns;
  for(int i = 0;i < iso_pair_count;i++){
    std::vector<std::vector<int> > temp_ns = umi_edit(umi,
                                            iso_index[iso_src(i,0)],
                                            iso_index[iso_src(i,1)],
                                            edit_thresh,k);

    if(temp_ns.size() > 0){
      umi_ns.insert(umi_ns.end(),temp_ns.begin(),temp_ns.end());
    }
  }
  return(umi_ns);
}
*/
