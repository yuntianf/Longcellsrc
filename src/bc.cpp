#include "bc.h"
#include "edit.h"

double mean(std::vector<double> num) {
  double n = num.size(), sum = 0;
  for (int i = 0; i < n; i++) {
    sum += num[i];
  }
  return(sum / n);
}

double var(std::vector<double> num) {
  double num_m = mean(num), sum = 0, n = num.size();
  for (int i = 0; i < n; i++) {
    sum += (num[i] - num_m) * (num[i] - num_m);
  }
  return(sum / (n - 1));
}

double update_sigma(std::vector<double> num, double sigma_start) {
  double n = num.size();
  if(n <= 1){
    return(sigma_start);
  }
  return(std::sqrt(var(num) * (n + 1) / (n + 0.5)) + sigma_start / n);
}

double update_prob(std::vector<double> m, double n) {
  int m_size = m.size();
  double sum = 0;
  for (int i = 0; i < m_size; i++) {
    sum += m[i];
  }
  return(sum / (n * m_size));
}

std::vector<int> kmer_include(std::string seq, std::set<std::string> dic) {
    std::vector<int> sc;
    int id = 0;
    for (std::set<std::string>::iterator it = dic.begin(); it != dic.end(); it++) {
        if (seq.find(*it) != seq.npos) {
            sc.push_back(id);
        }
        id++;
    }
    return(sc);
}

std::vector<std::vector<int> > barcodes_cos_vec(std::vector<std::string> barcodes, set<std::string> dic) {
    std::vector<std::vector<int> > bar_kmer_vec;
    int n = barcodes.size();

    for (int i = 0; i < n; i++) {
        int id = 0;
        std::vector<int> temp_vec = kmer_include(barcodes[i],dic);
        bar_kmer_vec.push_back(temp_vec);
    }
    return(bar_kmer_vec);
}

double cos_sim(std::vector<int> a, std::vector<int> b) {
    if (a.size() == 0 || b.size() == 0) {
        return(0);
    }
    std::vector<int> inter;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(inter));

    return(inter.size() / (sqrt(a.size()) * std::sqrt(b.size())));
}

bool cos_sim_comp(std::pair <int, double> cos1, std::pair <int, double> cos2) {
    return(cos1.second > cos2.second);
}

std::vector<std::string> barcode_cand_cos(std::string seq, std::vector<std::string> barcodes,
                                          std::set<std::string> dic, std::vector<std::vector<int> > index,
    int top, double thresh) {

    int bar_size = barcodes.size();

    std::vector<pair<int, double> > bar_cand;

    std::vector<int> sc = kmer_include(seq, dic);
    for (int j = 0; j < bar_size; j++) {
        double cos = cos_sim(sc,index[j]);
        if (cos >= thresh) {
            bar_cand.push_back(std::pair<int, double>(j, cos));
        }
    }
    std::sort(bar_cand.begin(), bar_cand.end(), cos_sim_comp);
    int count = bar_cand.size();
    int limit = std::min(top, count);

    std::vector<std::string> out;

    for (int i = 0; i < limit; i++) {
        out.push_back(barcodes[bar_cand[i].first]);
    }
    return(out);
}

std::vector<int> pos_filter(std::vector<double> start, std::vector<int> edit) {
    int num = start.size();

    std::vector<double> correct;
    std::vector<int> correct_id;

    int edit_count[16] = {}, offset_count[100] = {};
    int max_edit = 0, max_offset = 0;

    double edit_sum = 0;

    for (int i = 0; i < num; i++) {
        if (edit[i] == 0) {
            correct.push_back(start[i]);
            correct_id.push_back(i);
        }
        if (edit[i] > max_edit) {
            max_edit = edit[i];
        }
        edit_count[edit[i]]++;
        edit_sum += edit[i];
    }

    double edit_mean = edit_sum / num;

    if(correct.size() < 20){
      warning("Too few reads identified with confident cell barcode, please check if the barcode source is correct ");
      return(correct_id);
    }
    double mu = mean(correct);

    for (int i = 0; i < num; i++) {
        start[i] = std::round(abs(start[i] - mu));
        if (start[i] > max_offset) {
            max_offset = start[i];
        }
        offset_count[int(start[i])]++;
    }

    double edit_ratio[16] = {}, offset_ratio[100] = {};
    for (int i = max_edit; i >= 0; i--) {
        double count = 0;
        for (int j = i; j <= max_edit; j++) {
            count += edit_count[j];
        }
        edit_ratio[i] = count / double(num);
    }
    for (int i = max_offset; i >= 0; i--) {
        double count = 0;
        for (int j = i; j <= max_offset; j++) {
            count += offset_count[j];
        }
        offset_ratio[i] = count / double(num);
    }

    std::vector<std::pair<int, int> > preserve;
    for (int i = 0; i <= max_edit; i++) {
        for (int j = 0; j <= max_offset; j++) {
            if (std::pow(edit_ratio[i], edit_mean) * std::pow(offset_ratio[j], 0.5) >= 0.05 || i == 0) {
                preserve.push_back(std::pair<int, int>(i, j));
            }
        }
    }

    std::vector<int> preserve_id;
    for (int i = 0; i < preserve.size(); i++) {
        for (int j = 0; j < num; j++) {
            if (edit[j] == preserve[i].first && start[j] == preserve[i].second) {
                preserve_id.push_back(j);
            }
        }
    }
    std::sort(preserve_id.begin(), preserve_id.end());
    return(preserve_id);
}

// [[Rcpp::export]]
DataFrame barcodeMatch(std::vector<std::string> seq, std::vector<std::string> barcodes,
    double mu, double sigma, double sigma_start, int k, int batch,
    int top, double cos_thresh, double alpha, int edit_thresh,
    const int UMI_len, const int flank) {

    printf("start to build index\n");
    std::set<std::string> dic = kmer(barcodes, k,1);
    std::vector<std::vector<int> > index = barcodes_cos_vec(barcodes, dic);
    printf("index finished!\n");
    cout << "The size of index is " << dic.size() << endl;

    int seq_size = seq.size();
    int bar_size = barcodes.size();

    int bar_len = barcodes[1].length();

    int times = seq_size / batch;
    int seq_s, seq_e;

    std::vector<int> result_id;
    std::vector<std::string> result_bar;
    std::vector<double> result_pos;
    std::vector<int> result_edit;

    for (int i = 1; i <= times + 1; i++) {
      if(i % 100 == 0){
        checkUserInterrupt();
      }
        seq_s = batch * (i - 1);
        seq_e = batch * i;

        if (i == times + 1) {
            if (seq_size % batch) {
                seq_e = seq_size;
            }
            else {
                break;
            }
        }
        double interval_s = R::qnorm(alpha / 2, mu, sigma,true,false) < 0 ? 0 : (R::qnorm(alpha / 2, mu, sigma,true,false));;
        double interval_e = R::qnorm(1-alpha / 2, mu, sigma,true,false)+bar_len;

        interval_s = std::round(interval_s);
        interval_e = std::round(interval_e);

        for (int j = seq_s; j < seq_e; j++) {
            if (seq[j].length() < interval_s+bar_len-1) {
                continue;
            }
            std::string sub_seq = seq[j].substr(interval_s, interval_e - interval_s + 1);

            if (sub_seq.length() < bar_len) {
                continue;
            }
            else {
                std::vector<std::string> retain_barcodes = barcode_cand_cos(sub_seq, barcodes, dic, index, top, cos_thresh);
                int retain_barcodes_size = retain_barcodes.size();

                if (retain_barcodes_size == 0) {
                    continue;
                }
                int now = bar_len;
                pair<int, int> best(now, -100);
                int bar_id = 0;
                bool flag = 0;
                for (int p = 0; p < retain_barcodes_size; p++) {
                    std::pair<int, int> temp = minEditDist(sub_seq, retain_barcodes[p]);
                    if (temp.first < now) {
                        flag = 1;
                        bar_id = p;
                        best = temp;
                        now = temp.first;
                    }
                    else if (temp.first == now) {
                        if (abs(temp.second - mu) < abs(best.second - mu)) {
                            best = temp;
                            bar_id = p;
                        }
                    }
                }

                if (!flag) {
                    continue;
                }
                else if (best.first <= edit_thresh) {
                    result_id.push_back(j);
                    result_bar.push_back(retain_barcodes[bar_id]);
                    result_pos.push_back(interval_s + best.second);
                    result_edit.push_back(best.first);
                }
            }

        }
        if(result_pos.size() > batch && result_pos.size() < 10000){
          mu = mean(result_pos);
          sigma = update_sigma(result_pos, sigma_start);
        }
    }

    std::vector<int> out_id;
    std::vector<std::string> out_bar;
    std::vector<int> out_pos;
    std::vector<int> out_edit;

    //cout << result_e.size() << endl;
    //DataFrame out = DataFrame::create(Named("end") =  result_e, Named("edit") = result_edit);
    //return(out);
    if(result_pos.size() > 0){
      //std::vector<int> preserve_id = pos_filter(result_e, result_edit);
      std::vector<int> preserve_id;
      for(int i = 0;i < result_pos.size();i++){
        preserve_id.push_back(i);
      }
      int preserve_num = preserve_id.size();
      cout << "There are " << preserve_num << " sequences identified with a barcode" << endl;

      for (int i = 0; i < preserve_num; i++){
        out_id.push_back(result_id[preserve_id[i]]);
        out_bar.push_back(result_bar[preserve_id[i]]);
        out_pos.push_back(int(result_pos[preserve_id[i]]));
        out_edit.push_back(result_edit[preserve_id[i]]);
      }

    }
    std::vector<std::vector<std::string> > neighbor = NeighborExtract(seq,out_id,out_pos,
                                                                      UMI_len,flank,bar_len);
    DataFrame out = DataFrame::create(Named("id") =  out_id, Named("barcode") = out_bar,
                                      Named("pos") = out_pos,Named("edit") = out_edit,
                                      Named("umi") = neighbor[0],Named("adapter") = neighbor[1]);

    //List out = List::create(Named("mu") = Mu, Named("var") = Var,Named("count") = cand);
    return(out);
}

// [[Rcpp::export]]
std::vector<std::vector<std::string> > NeighborExtract(std::vector<std::string> reads,std::vector<int> id,
                                                       std::vector<int> pos,const int UMI_len,
                                                       const int flank,const int bar_len){
  int len = id.size();
  try{
    if(len != pos.size()){
      throw "The size of reads is and the start position don't match!";
    }
  }
  catch (const char* errorMsg) {
    std::cerr << "Error: " << errorMsg << std::endl;
  }
  //std::vector<int> id_vector;
  std::vector<std::string> UMI_vector;
  std::vector<std::string> adapter_vector;
  std::string umi = "";
  std::string adapter = "";
  for(int i = 0;i < len;i++){
    umi = "";
    adapter = "";
    if(pos[i] >= UMI_len+flank){
      umi = reads[id[i]].substr(pos[i]-UMI_len-flank,UMI_len+2*flank);
    }
    if(pos[i]+bar_len+UMI_len+2*flank <= reads[id[i]].size()){
      adapter = reads[id[i]].substr(pos[i]+bar_len,UMI_len+2*flank);
    }
    UMI_vector.push_back(umi);
    adapter_vector.push_back(adapter);

  }

  std::vector<std::vector<std::string> > out;
  out.push_back(UMI_vector);
  out.push_back(adapter_vector);
  return(out);
}

