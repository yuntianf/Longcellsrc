#include "tag_extraction.h"

/*
// [[Rcpp::export]]
std::vector<std::string> extractSoftclip(std::string seq,StringVector mark,NumericVector count,
                                         const std::string strand,int toolkit,
                                         const int search_len,const int polyA_bin,
                                         const int polyA_count,const int polyA_len){
  // extract softclips at both ends
  std::string softclip_5 = "";
  std::string softclip_3 = "";

  int len = mark.size();
  if(len > 0){
    if(mark[0] == "S"){
      int start = 0;
      int end = count[0];
      softclip_5 = seq.substr(start,end-start);
    }
    if(mark[len-1] == "S"){
      int start = seq.size() - count[len-1] + 1;
      softclip_3 = seq.substr(start);
    }
  }

  // extract softclips at both ends
  std::string search_seq = "";
  bool polyA_exist = false;
  if(softclip_5.size() > 0 || softclip_3.size() > 0){
    int site = 0;
    if(strand == "+"){
      polyA_exist = polyADetect(softclip_3,polyA_bin,polyA_count);
      if(toolkit == 5){
        site = std::max(0,int(softclip_5.size()-search_len + 1));
        search_seq = softclip_5.substr(site);
      }
      else if(toolkit == 3){
        softclip_3 = reverseComplement(polyARm(softclip_3,polyA_len));
        site = std::max(0,int(softclip_3.size()-search_len + 1));
        search_seq = softclip_3.substr(site);
      }
    }
    else{
      polyA_exist = polyADetect(reverseComplement(softclip_5),polyA_bin,polyA_count);
      if(toolkit == 5){
        softclip_3 = reverseComplement(softclip_3);
        site = std::max(0,int(softclip_3.size()-search_len + 1));
        search_seq = softclip_3.substr(site);
      }
      else if(toolkit == 3){
        softclip_5 = polyARm(reverseComplement(softclip_5),polyA_len);
        softclip_5 = reverseComplement(softclip_5);
        site = std::max(0,int(softclip_5.size()-search_len + 1));
        search_seq = softclip_5.substr(site);
      }
    }
  }
  std::vector<std::string> out = {search_seq,std::to_string(polyA_exist)};
  return(out);
}
 */

std::string filePrefix(const std::string& path) {
  // Find the last occurrence of '/' or '\'
  size_t last_slash = path.find_last_of("/\\");
  if (last_slash == std::string::npos) {
    last_slash = -1; // No directory separator found
  }
  // Find the last occurrence of '.'
  size_t last_dot = path.find_last_of('.');
  if (last_dot == std::string::npos) {
    last_dot = path.length(); // No extension found
  }
  // Extract file name prefix (name without extension)
  return path.substr(last_slash + 1, last_dot - last_slash - 1);
}

void writeSeq(vector<KSeq> records,string path){
  SeqStreamOut oss(path.c_str(),true, format::fastq);

  for(auto record:records){
    oss << record;
  }
}

// [[Rcpp::export]]
void fastqSplit(const char* fastq_path,const char* out_path,const int batch){

  KSeq record;
  SeqStreamIn iss(fastq_path);

  std::string input(fastq_path);
  std::string prefix = filePrefix(input);

  std::string output(out_path);
  std::string suffix = ".fq.gz";
  std::string sep = "/";
  std::string split = "_";

  int id = 0,i = 0;
  std::string out_name = "";
  vector<KSeq> container;

  while (iss >> record){
    container.push_back(record);
    i++;
    if(i >= batch){
      i = 0;
      out_name = (output+sep+prefix+split+to_string(id)+suffix);
      writeSeq(container,out_name);
      container.clear();
      id++;
    }
  }
}

//' extractTagFastq
 //'
 //' This function extract the tag region for a fastq.
 //'
 //' @inheritParams extractTagFastq
 //' @param fastq_path The filename of the input fastq.
 //' @param out_path The filename for the output polished fastq.
 //' @return A dataframe including the read name, the tag region and the polyA existence for each read.
 //' @export
// [[Rcpp::export]]
DataFrame extractTagFastq(const char* fastq_path,const char* out_path,
                                         const std::string adapter, const int toolkit,
                                         const int window, const int step,
                                         const int left_flank, const int right_flank, const bool drop,
                                         const int polyA_bin,
                                         const int polyA_base_count,const int polyA_len){

  KSeq record;
  SeqStreamIn iss(fastq_path);
  SeqStreamOut oss(out_path,true, format::fastq);

  std::vector<std::string> names;
  std::vector<std::string> tags;
  std::vector<std::string> polyAs;

  std::vector<std::string> tag;
  int i = 0, j = 0;
  while (iss >> record) {
    i++;
    if(i % 100 == 0){
      checkUserInterrupt();
    }
    tag = extractTag(record,adapter,toolkit,window,step,
                     left_flank, right_flank, drop,
                     polyA_bin, polyA_base_count, polyA_len);
    if(tag[4] == "1"){
      j++;
    }
    //cout << tag[0] << endl << tag[1] << endl << tag[2] << endl << tag[3] << endl;
    if(tag[0].size() >= 30){
      names.push_back(record.name);
      tags.push_back(tag[0]);

      record.seq = tag[1];
      record.qual = tag[2];

      polyAs.push_back(tag[3]);

      oss << record;
    }
  }
  Rcpp::Rcout << "In total " << j << " reads out of " << i << " reads are identified with a valid adapter." << endl;
  Rcpp::Rcout << "Within them, " << tags.size() << " reads are identified with a valid tag region." << endl;
  DataFrame out = DataFrame::create(Named("name") = names, Named("tag") = tags, Named("polyA") = polyAs);
  return(out);
}

//' extractTag
 //'
 //' This function extract the tag region for a read
 //'
 //' @inheritParams extractTag
 //' @param polyA_bin The window to search for polyA.
 //' @param polyA_base_count The minimum threshold for the times that the base appears in the sequence.
 //' @param polyA_len The maximum length of polyA to be preserved in the trimmed cDNA.
 //' @param flank The size of the flank sequence when extract the UMI
 //' @return A string vector including read name, trimmed sequence, trimmed read quality and
 //' polyA existence for a read.
std::vector<std::string> extractTag(KSeq record, const std::string adapter,
                                    const int toolkit,
                                    const int window, const int step,
                                    const int left_flank, const int right_flank, const bool drop,
                                    const int polyA_bin, const int polyA_base_count,const int polyA_len){
  std::string seq = record.seq;
  std::string qual = record.qual;

  //cout << seq.size() << " " << qual.size() << endl;
  int pos = -1, rpos = -1;
  std::string rseq = reverseComplement(seq);
  if(toolkit == 5){
    pos = strSlideSearch(seq,adapter,window,step,false);
    rpos = strSlideSearch(rseq,adapter,window,step,false);
  }
  else if(toolkit == 3){
    pos = strSlideSearch(seq,adapter,window,step,true);
    rpos = strSlideSearch(rseq,adapter,window,step,true);
  }

  //std::cout << pos << ":" << rpos << endl;

  std::string rqual(qual);
  std::reverse(rqual.begin(), rqual.end());

  std::string tag = "";
  bool flag = 1;
  bool polyA = 0;

  if(pos == -1 && rpos == -1){
    flag = 0;
  }
  else{
    if(pos == -1 && rpos != -1){
      seq = rseq;
      pos = rpos;
      qual = rqual;
    }

    const char base = 'A';
    bool polyA = polyADetect(seq,polyA_bin,polyA_base_count,base);

    if(pos != -1){
      if(drop){
        string left = seq.substr(std::max(0,pos-left_flank),min(pos,left_flank));
        string right = seq.substr(std::min(pos+adapter.size(),seq.size()),right_flank);
        tag = left+right;
      }
      else{
        tag = seq.substr(std::max(0,pos-left_flank),std::min(pos,left_flank)+adapter.size()+right_flank);
      }
      if(toolkit == 5){
        int start = std::min(seq.size(),pos+adapter.size()+right_flank);
        seq = seq.substr(start);
        qual = qual.substr(start);
      }
      else if(toolkit == 3){
        int end = std::max(0,pos-left_flank);
        seq = seq.substr(0,end);
        qual = qual.substr(0,end);
      }
    }
  }



  /*
  if(pos != -1 && rpos == -1){
    if(toolkit == 5){
      tag = seq.substr(std::max(0,pos-len),min(pos+2,len+2));
      seq = seq.substr(pos);
      qual = qual.substr(pos);
    }
    else if(toolkit == 3){
      tag = seq.substr(std::max(0,pos-len),min(pos+2,len+2));
      seq = seq.substr(0,std::max(0,pos-30));
      qual = qual.substr(0,std::max(0,pos-30));
    }

  }
  else if(pos == -1 && rpos != -1){
    if(toolkit == 5){
      tag = rseq.substr(std::max(0,rpos-len),min(rpos+2,len+2));
      seq = rseq.substr(rpos);
      qual = rqual.substr(rpos);
    }
    else if(toolkit == 3){
      tag = rseq.substr(std::max(0,rpos-len),min(rpos+2,len+2));
      seq = rseq.substr(0,std::max(0,rpos-30));
      qual = rqual.substr(0,std::max(0,rpos-30));
    }
  }
  */

  //cout << seq.size() << " " << qual.size() << endl;

  std::vector<std::string> result;

  //cout << seq.size() << " " << qual.size() << endl;

  result.push_back(tag);
  result.push_back(seq);
  result.push_back(qual);
  result.push_back(std::to_string(polyA));
  result.push_back(std::to_string(flag));

  return(result);
}



//' strSlideSearch
 //'
 //' This function is to search for a substring in a string in a slide window way
 //'
 //' @inheritParams strSubset
 //' @param record A KSeq object recording a read in the fastq file.
 //' @param adapter A string to be searched in the read.
 //' @return An int to indicate the position of the adapter in the sequence, -1 if not found.
 //' @export
 // [[Rcpp::export]]
int strSlideSearch(std::string seq,const std::string adapter,
                   const int window, const int step,const bool first){
  int pos = -1;
  if(first){
    pos = seq.find(adapter);
  }
  else{
    pos = seq.rfind(adapter);
  }

  // cout << pos << endl;
  if(pos != string::npos){
    return(pos);
  }
  else{
    std::vector<int> pos_vec;
    std::vector<std::string> sub_adapter = strSubset(adapter,window,step);
    for(auto i:sub_adapter){
      if(first){
        pos = seq.find(i);
      }
      else{
        pos = seq.rfind(i);
      }
      if(pos != string::npos){
        pos_vec.push_back(pos);
      }
    }
    if(pos_vec.size() == 0){
      return(-1);
    }
    else{
      if(pos_vec.size() < sub_adapter.size()/2){
        return(-1);
      }
      const int n = pos_vec.size();
      return(std::accumulate(pos_vec.begin(), pos_vec.end(),0) / n);
    }
  }
}

//' strSubset
 //'
 //' This function to subset a string in a slide window way
 //'
 //' @inheritParams strSubset
 //' @param str A string to be subset.
 //' @param window The size of the substring.
 //' @param step The distance between two substrings.
 //' @return A string vector inlcuding all substrings.
std::vector<std::string> strSubset(std::string str,const int window, const int step){
  std::vector<std::string> vec;
  int len = str.size();
  if(window >= len){
    vec.push_back(str);
  }
  else{
    int start = 0;
    int end = window;
    while(end < len){
      std::string sub = str.substr(start,window);
      vec.push_back(sub);
      start += step;
      end += step;
    }
    if(end != len){
      std::string sub = str.substr(len-window+1,window);
      vec.push_back(sub);
    }
  }
  return(vec);
}

//' polyADetect
 //'
 //' This function detect if polyA exists in the cDNA.
 //'
 //' @param seq A DNA seq string.
 //' @param bin The window to search for polyA.
 //' @param count The minimum threshold for the times that the base appears in the sequence.
 //' @param base The base to be searched, should be 'A' for polyA detection.
 //' @return A bool to indicate if the polyA exist.
bool polyADetect(std::string seq,const int bin, const int count,const char base){
  int len = seq.size();

  int i = 0;
  while(i <= (len-bin)){
    std::string subseq = seq.substr(i,bin);
    int bc = baseCount(subseq, base);

    if(bc < count){
      int next = count - bc;
      i = i + next;
    }
    else{
      return(true);
    }
  }

  return(false);
}


//' polyARm
 //'
 //' This function trim the polyA in a cDNA sequence.
 //'
 //' @param seq A DNA seq string.
 //' @param polyA_len The maximum length of polyA to be preserved in the trimmed cDNA.
 //' @return A string as the trimmed cDNA.
size_t polyARm(std::string seq, const int polyA_len){
  std::string polyA = replicate("A",polyA_len);
  size_t pos = seq.rfind(polyA);

  //if (pos != std::string::npos){
  //  seq = seq.substr(pos+polyA_len-5);
  //}
  return(pos);
}


//' baseCount
 //'
 //' This function count the time that a base appearing in a DNA seq.
 //'
 //' @param seq A DNA seq string.
 //' @param base A character.
 //' @return A int to indicate the time that this base appears.
// [[Rcpp::export]]
int baseCount(std::string seq, char base){
  int len = seq.size();
  int count = 0;
  for(int i = 0;i < len;i++){
    if(seq[i] == base){
      count += 1;
    }
  }
  return(count);
}

//' reverseComplement
 //'
 //' This function generate the reverse complement version for a DNA seq
 //'
 //' @param sequence A DNA seq string.
 //' @return A string as reverse complemented DNA.
std::string reverseComplement(const std::string& sequence) {
  std::string complement;
  complement.reserve(sequence.length());

  for (char base : sequence) {
    switch (base) {
    case 'A':
      complement += 'T';
      break;
    case 'T':
      complement += 'A';
      break;
    case 'G':
      complement += 'C';
      break;
    case 'C':
      complement += 'G';
      break;
    default:
      complement += base;
    }
  }

  std::reverse(complement.begin(), complement.end());
  return(complement);
}

//' replicate
 //'
 //' This function generate strings with replicated patterns
 //'
 //' @param mode A string to be replicated.
 //' @param times An int to indicate the number of replicates.
 //' @return A string with replicated pattern.
std::string replicate(std::string mode,int times){
  std::string out = "";
  for(int i = 0;i < times;i++){
    out = out+mode;
  }
  return(out);
}


