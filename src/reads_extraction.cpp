#include "reads_extraction.h"

bool CHOP = false;

// [[Rcpp::export]]
DataFrame extractReads(std::vector<std::string> seq,std::vector<std::string> cigar,
                       NumericVector pos,NumericMatrix annotation,const std::string strand,
                       int toolkit,int end_flank,int splice_site_bin){
  std::vector<int> id;
  std::vector<std::string> isoformVec;
  std::vector<int> isoEndVec;
  std::vector<std::string> polyA;

  int refer_start = annotation(0,0);
  int refer_end = annotation(annotation.nrow()-1,1);

  int len = seq.size();
  std::string read_polyA = "1";
  for(int i =0;i < len;i++){
    List cigar_list = cigarProcess(cigar[i]);
    StringVector mark = cigar_list[0];
    NumericVector count = cigar_list[1];
    if(!cigarcheck(mark,count)){
      continue;
    }

    //int end = seqEnd(pos[i],mark,count);

    NumericMatrix isoform = extractIsoform(pos[i],mark,count,refer_start,refer_end,end_flank);
    isoform = splicesiteCorrect(isoform,annotation,strand,splice_site_bin);
    if(isoform(0,0) == -1){
      continue;
    }
    read_polyA = "1";
    if(CHOP){
      read_polyA = "0";
      CHOP = false;
    }
    int isoEnd= isoformEnd(isoform,strand);
    std::string isoformSeq = isoform2string(isoform);


    if(isoformSeq.size() > 0){
      id.push_back(i+1);
      isoformVec.push_back(isoformSeq);
      isoEndVec.push_back(isoEnd);
      polyA.push_back(read_polyA);
    }
  }

  DataFrame out = DataFrame::create(Named("id") =  id,
                                    Named("isoform") = isoformVec,
                                    Named("isoend") = isoEndVec,
                                    Named("polyA") = polyA);
  return(out);
}

NumericMatrix extractIsoform(int read_start,
                             StringVector mark, NumericVector count,
                             int refer_start, int refer_end,
                             int flank){
  NumericVector blocks;

  int start = read_start;
  int end = start;
  int len = mark.size();

  bool flag = false;
  for(int i = 0;i < len;i++){
    if(mark[i] == "M" || mark[i] == "D"){
      end += count[i];
    }
    else if(mark[i] == "N"){
      if(start >= refer_start-flank && end <= refer_end+flank){
        blocks.push_back(start);
        blocks.push_back(end-1);
      }
      else{
        flag = true;
      }
      start = end + count[i];
      end = start;
    }
  }
  if(start >= refer_start-flank && end <= refer_end+flank){
    blocks.push_back(start);
    blocks.push_back(end-1);
  }
  else{
    flag = true;
  }
  if(!flag){
    NumericMatrix isoform(2,blocks.size()/2,blocks.begin());
    return(transpose(isoform));
  }
  else{
    NumericVector temp = NumericVector::create(-1, -1);
    NumericMatrix isoform(1, 2, temp.begin());
    return(isoform);
  }
}

NumericMatrix chopOutBound(NumericMatrix isoform,NumericMatrix annotation,string strand){
  if(isoform(0,0) == -1){
    return(isoform);
  }

  int start = 0, end = isoform.nrow()-1;

  while(isoform(start,1) < annotation(0,0)){
    start++;
  }
  while(isoform(end,0) > annotation(annotation.nrow()-1,1)){
    end--;
  }

  if(strand == "+" && end < isoform.nrow()-1){
    CHOP = true;
  }
  if(strand == "-" && start > 0){
    CHOP = true;
  }
  if(start <= end){
    return(isoform(Range(start,end),_));
  }
  else{
    NumericVector temp = NumericVector::create(-1, -1);
    NumericMatrix out(1, 2, temp.begin());
    return(out);
  }
}

NumericMatrix splicesiteCorrect(NumericMatrix isoform,NumericMatrix annotation,
                                string strand,int bin){
  if(isoform.nrow() == 0){
    return(isoform);
  }
  isoform = chopOutBound(isoform,annotation,strand);

  if(isoform(0,0) == -1){
    return(isoform);
  }
  if(isoform.nrow() == 1){
    if(isoform(0,0) < annotation(0,0)){
      isoform(0,0) = annotation(0,0);
    }
    if(isoform(0,1) > annotation(annotation.nrow()-1,1)){
      isoform(0,1) = annotation(annotation.nrow()-1,1);
    }
    return(isoform);
  }

  NumericVector exon_start = isoform(_,0);
  NumericVector exon_end = isoform(_,1);

  NumericVector refer_start = annotation(_,0);
  NumericVector refer_end = annotation(_,1);

  int i_len = isoform.nrow();
  int r_len = annotation.nrow();

  if(exon_start(0) < refer_start(0)){
    exon_start(0) = refer_start(0);
  }
  if(exon_end(i_len-1) > refer_end(r_len-1)){
    exon_end(i_len-1) = refer_end(r_len-1);
  }

  int i = 1,j = 0;
  while(i < i_len){
    while(exon_start[i] > refer_start[j]+bin){
      if(j == r_len){
        break;
      }
      j++;
    }
    int mindis = bin*10;
    int cand = exon_start[i];
    while(exon_start[i] <= refer_start[j]+bin && exon_start[i] >= refer_start[j]-bin){
      int dis = abs(exon_start[i] - refer_start[j]);
      if(dis < mindis){
        mindis = dis;
        cand = refer_start[j];
      }
      j++;
    }
    exon_start[i] = cand;
    i++;
  }

  i = 0,j = 0;
  while(i < (i_len-1)){
    while(exon_end[i] > refer_end[j]+bin){
      if(j == r_len-1){
        break;
      }
      j++;
    }
    int mindis = bin*10;
    int cand = exon_end[i];
    while(exon_end[i] <= refer_end[j]+bin && exon_end[i] >= refer_end[j]-bin){
      int dis = abs(exon_end[i] - refer_end[j]);
      if(dis < mindis){
        mindis = dis;
        cand = refer_end[j];
      }
      j++;
    }
    exon_end[i] = cand;
    i++;
  }

  NumericMatrix out(i_len,2);
  out(_,0) = exon_start;
  out(_,1) = exon_end;

  return(out);
}

int isoformEnd(NumericMatrix isoform, std::string strand){
  int len = isoform.nrow();
  if(strand == "+"){
    return(int(isoform(len-1,1)));
  }
  else{
    return(int(isoform(0,0)));
  }
}

std::string isoform2string(NumericMatrix isoform,std::string sep){
  int len = isoform.nrow();

  std::string seq = "";
  for(int i = 0;i < len;i++){
    seq = seq + std::to_string(int(isoform(i,0))) + "," + std::to_string(int(isoform(i,1)));
    if(i != len-1){
      seq += sep;
    }
  }
  return(seq);
}

// [[Rcpp::export]]
List cigarProcess(std::string cigar){
  std::regex pattern(R"(([0-9]+(?:\.[0-9]+)?)\s*([a-zA-Z]+))");

  std::sregex_iterator iter(cigar.begin(), cigar.end(), pattern);
  std::sregex_iterator end;

  StringVector mark;
  NumericVector count;

  for (; iter != end; ++iter) {
    std::smatch match = *iter;

    std::string numberstring = match[1].str();
    std::string cigarstring = match[2].str();

    int number = stoi(numberstring);

    mark.push_back(cigarstring);
    count.push_back(number);
  }

  if(mark.size() != count.size()){
    printf("Warning: the size of the mark and count from cigar don't match!");
  }
  List out = List::create(mark, count);
  return(out);
}

bool cigarcheck(StringVector mark,NumericVector count){
  int size = count.size();
  for(int i = 0;i < size;i++){
    if(mark[i] == "I" || mark[i] == "D"){
      if(count[i] > 6){
        return(false);
      }
    }
  }
  return(true);
}

bool within(int posA,int posB,const int bin){
  return(posA >= posB-bin && posA <= posB+bin);
}
// [[Rcpp::export]]
int seqEnd(int start,StringVector mark,NumericVector count){
  int len = mark.size();
  int end = start;
  for(int i = 0;i < len;i++){
    if(mark[i] == "M" || mark[i] == "D" || mark[i] == "N"){
      end += count[i];
    }
  }
  return(end);
}
