/*
 * bcp_algo.cpp
 *
 *  Created on: Sep 27, 2014
 *      Author: Xin Feng 
 */
#include "bcp_algo.h"
#include "common/boost_header.h"
#include "common/stl_header.h"
#include "utils/assert_helpers.h"
#include "utils/logger.h"
#include "utils/debug.h"
#include "short_reads/readstools.h"
#include "region_detector/calledpeak.h"
#include "utils/Tracer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "cppoisson_HM.h"


namespace {

  int cmp(const void *a, const void *b) {
    return *(int *) a - *(int *) b;
  }

}

using namespace std;
using namespace TNT;
using namespace boost;

void bcp_algo::loadPosData(data_t& data, Reads& reads) {

  //Assume reads in Reads already normalized by chrs in the main app.cpp
  foreach(std::string chr, reads.pos_reads.chrs()) {
    std::vector<uint32_t>::iterator it = reads.pos_reads.begin_of(chr);
    while (it != reads.pos_reads.end_of(chr)) {
      // This is safe due to the length of chromosome
      data[chr2Num[chr]].push_back((int) (*it++));
    }
  }
}

void bcp_algo::loadNegData(data_t& data, Reads& reads) {

  //Assume reads in Reads already normalized by chrs in the main app.cpp
  foreach(string chr, reads.neg_reads.chrs()) {
    std::vector<uint32_t>::iterator it = reads.neg_reads.begin_of(chr);
    while (it != reads.neg_reads.end_of(chr)) {
      // This is safe due to the length of chromosome
      data[chr2Num[chr]].push_back((int) (*it++));
    }
  }
}

void bcp_algo::insertPeak(const string& chr, called_peak& pk) {
    map<string, vector<called_peak> >::iterator it;
    it = _resultRegions.find(chr);

    if (it == _resultRegions.end()) {
        vector<called_peak> peaks;
        _resultRegions[chr] = peaks;
    }
    _resultRegions[chr].push_back(pk);
}

void bcp_algo::cmain(Reads& treads, Reads& creads, cmd_option_parser& option) {
  utils::TimeStampTracer tracer(std::cout, option.getVerboseRequested());

  p_value = option.getCut_off();
  win_size = option.slidingWinSize;
  frag_size = option.getExt_length();

  data_t Plus_data(N), Minus_data(N);
  data_t Plus_input(N), Minus_input(N);

  int *temp;
  int pre,num,m,l1,l2,t1,t2,seg_len;
  double num_temp, thre, c, max, max1, lambda;
  num_win = 0;


  buildChrMap(treads);
  loadPosData(Plus_data, treads);
  loadNegData(Minus_data, treads);
  loadPosData(Plus_input, creads);
  loadNegData(Minus_input, creads);

  std::pair<uint32_t, std::string> mi;
  uint32_t i = 0;
  foreach (mi, num2Chr) {
    i = mi.first;
    /*The data part*/
    data_frag = new int[Plus_data[i].size() + Minus_data[i].size() + 1];
    num = 1;
    data_frag[0] = 0;
    if (Plus_data[i].size() > 0) {
      temp = new int[Plus_data[i].size()];

      for (unsigned int j = 0; j < Plus_data[i].size(); j++)
        temp[j] = Plus_data[i][j];

      qsort(temp, Plus_data[i].size(), sizeof(int), cmp);

      pre = -1;
      for (unsigned int j = 0; j < Plus_data[i].size(); j++) {
        if (temp[j] != pre) {
          data_frag[num] = temp[j];
          Add(temp[j], temp[j] + frag_size - 1);
          num++;
        }
        pre = temp[j];
      }
      delete[] temp;
    }

    if (Minus_data[i].size() > 0) {
      temp = new int[Minus_data[i].size()];

      for (unsigned int j = 0; j < Minus_data[i].size(); j++)
        temp[j] = Minus_data[i][j];

      qsort(temp, Minus_data[i].size(), sizeof(int), cmp);

      pre = -1;
      for (unsigned int j = 0; j < Minus_data[i].size(); j++) {
        if (temp[j] != pre) {
          data_frag[num] = temp[j];
          Add(temp[j], temp[j] + frag_size - 1);
          num++;
        }
        pre = temp[j];
      }
      delete[] temp;
    }
    // in fact the num_data should be num-1, just for convenient for qsort;
    num_data_frag = num;
    qsort(data_frag, num_data_frag, sizeof(int), cmp);

    //The input part
    //TODO: Refactor to extract the function
    input_frag = new int[Plus_input[i].size() + Minus_input[i].size() + 1];
    num = 1;
    input_frag[0] = 0;
    if (Plus_input[i].size() > 0) {
      temp = new int[Plus_input[i].size()];

      for (unsigned int j = 0; j < Plus_input[i].size(); j++)
        temp[j] = Plus_input[i][j];

      qsort(temp, Plus_input[i].size(), sizeof(int), cmp);

      pre = -1;
      for (unsigned int j = 0; j < Plus_input[i].size(); j++) {
        if (temp[j] != pre) {
          input_frag[num] = temp[j];
          num++;
        }
        pre = temp[j];
      }
      delete[] temp;
    }

    if (Minus_input[i].size() > 0) {
      temp = new int[Minus_input[i].size()];

      for (unsigned int j = 0; j < Minus_input[i].size(); j++)
        temp[j] = Minus_input[i][j];

      qsort(temp, Minus_input[i].size(), sizeof(int), cmp);

      pre = -1;
      for (unsigned int j = 0; j < Minus_input[i].size(); j++) {
        if (temp[j] != pre) {
          input_frag[num] = temp[j];
          num++;
        }
        pre = temp[j];
      }
      delete[] temp;
    }
    num_input_frag = num; 
    // in fact the num_input should be num-1, just for convenient for qsort;
    qsort(input_frag, num_input_frag, sizeof(int), cmp);
    if (Plus_data[i].size() + Minus_data[i].size() > 0) {
      win_data = new double[L2][4];
      trans2window(i);

      len_bf = sum_bf = 0.0;
      for (m = 2; m <= num_win; m++) {
        len_bf += win_data[m][3];
        sum_bf += win_data[m][2] * win_data[m][3];
      }
      average_bf = sum_bf / len_bf;

      Matrix<double> obs(num_win + 1, 4);
      Matrix<double> data(num_win + 1, 5);

      m = 0;
      obs[m][0] = obs[m][1] = obs[m][2] = obs[m][3] = 0.0;
      for (m = 1; m <= num_win; m++) {
        obs[m][0] = win_data[m][0];
        obs[m][1] = win_data[m][1];
        obs[m][2] = floor(win_data[m][2] + 0.5);
        obs[m][3] = win_data[m][3];
      }

      len_aft = sum_aft = 0.0;
      for (m = 2; m <= num_win; m++) {
        len_aft += obs[m][3];
        sum_aft += obs[m][2] * obs[m][3];
      }
      average_aft = sum_aft / len_aft;

      for (m = 1; m < num_win; m++)
        if ((obs[m][2] == 0) && (obs[m][3] == 1)
            && ((obs[m - 1][2] != 0) || (obs[m + 1][2] != 0)))
          obs[m][3] = 0;

      num_temp = floor(log10((double) num_win) + 0.5);
      double p = 1.0 / pow(10.0, num_temp);
      if (average_aft < 0.1 && p < 0.0001)
        p = p * 100;
      if (average_aft >= 0.1 && average_aft <= 0.5 && p < 0.0001)
        p = p * 10;

      if (win_size >= 200) {
        t1 = 10;
        t2 = 5;
      }
      if (win_size < 200 && win_size >= 100) {
        t1 = 7;
        t2 = 4;
      }
      if (win_size < 100) {
        t1 = 5;
        t2 = 3;
      }
      if (t1 >= num_win) {
        t1 = 5;
        t2 = 3;
      }
      if (t1 >= num_win) {
        tracer << "Too few data for the model. Program stopped prematurely.\n";
        exit(1);
      }

      cppoisson tmp(obs, p, 1.0, 1.0, t1, t2);
      tmp.BcmixSmooth();

      m = 0;
      data[m][0] = data[m][1] = data[m][2] = data[m][3] = data[m][4] =
        0.0;
      for (m = 1; m <= num_win; m++) {
        data[m][0] = obs[m][0];
        data[m][1] = obs[m][1];
        data[m][2] = obs[m][2];
        data[m][3] = obs[m][3];
        data[m][4] = tmp.estPara[m];
      }


      thre = 0.9;
      //calculationg the factorial first
      st[0] = 0;
      for (m = 1; m <= 100000; m++)
        st[m] = st[m - 1] + log(m);

      m = 0;
      while (prob_pois(m, average_aft) < thre)
        m++;

      if (m == 0)
        cutline = 1.0;
      if (m > 3)
        cutline = m * 1.0 - 0.5;
      if (m >= 1 && m <= 3)
        cutline = m * 1.0;

      ss = new double[num_win][6];

      seg(data, cutline);

      rseg = new double[num_seg + 1][7];

      int m1 = 0;
      for (m = 1; m <= num_allseg; m++) {
        if (ss[m][3] >= cutline) {
          m1++;
          rseg[m1][0] = ss[m][0];
          rseg[m1][1] = ss[m][1];
          rseg[m1][2] = ss[m][2];
          rseg[m1][3] = ss[m][3];
          rseg[m1][4] = rseg[m1][5] = 0;
          rseg[m1][6] = 1.0;
        }
      }

      c = (double) num_data_frag / (double) num_input_frag;

      max = (num_input_frag - 1)
        / (double) (input_frag[num_input_frag - 1] + frag_size - 1
            - input_frag[1]);

      for (m = 1; m <= num_seg; m++) {
        seg_len = (int) rseg[m][2];
        max1 = max * seg_len;

        l1 = frag_count1(m);
        rseg[m][4] = (double) l1;

        l2 = frag_count2(m);
        lambda = l2 >= max1 ? l2 : max1;

        rseg[m][5] = c * lambda;
        rseg[m][6] = fabs(1.0 - prob_pois(l1, rseg[m][5]));
      }

      for (m = 1; m <= num_seg; m++){
        if ((rseg[m][6] <= p_value) && (rseg[m][4] > rseg[m][5])){
        std::vector<uint32_t> tmp;
        called_peak pk((rseg[m][0] - 1), rseg[m][1], rseg[m][6] ,
        		rseg[m][6] , rseg[m][2], rseg[m][3], tmp);
        insertPeak(num2Chr[i], pk);
        }
      }

      tracer <<"Discovered "<<_resultRegions[num2Chr[i]].size()<<"\tregions in "<<num2Chr[i]<<".\n";

      delete[] rseg;
      delete[] ss;
      delete[] win_data;
    } else{
      tracer << num2Chr[i] << "is empty.\n" ;
    }

      for (m = 0; m < L1; m++)
        Weight[m] = 0;
      delete[] data_frag;
      delete[] input_frag;
  }
}

void bcp_algo::printhelp() {
  cout << "\t-1\tThe ChIP-seq data set you want to input." << endl;
  cout << "\t-2\tThe control/input data set you want to input." << endl;
  cout
    << "\t-f\tThe fragment size to which we extend the reads in pre-processing data step. Default:200bp."
    << endl;
  cout
    << "\t-w\tThe window size we apply the adjcaent window in pre-processing data step. Default:200bp."
    << endl;
  cout
    << "\t-p\tThe p_value you want to use for remove false positive based on control data. Range:1e-2---1e-6, default is 1e-3."
    << endl;
  cout << "\t-3\tThe results data set with 5 columns you want to output."
    << endl;
  cout
    << "\t-h\tThe flag indicates whether you want to output the help manual. 1 means displaying the help instructions while 0 means NOT. Default: 0. If you set the flag as 1, the program would not works except displaying the help manual."
    << endl;
}

void bcp_algo::Add(int aa, int bb) {
  aa++;
  int head_index = aa / win_size;
  int tail_index = bb / win_size;

  if (head_index == tail_index)
    Weight[head_index] += bb - aa + 1;
  else //if(head_index<tail_index)
  {

    Weight[head_index] += win_size - aa % win_size;

    for (int k = 1; k < tail_index - head_index; k++)
      Weight[head_index + k] += win_size;

    Weight[tail_index] += bb % win_size + 1;

  }
}

void bcp_algo::trans2window(int r) {
  int pre_index = -1;
  int i1 = 0;

  for (unsigned int j = 0; j < L1; j++) {
    if (Weight[j] > 0) {
      int n1 = j - pre_index - 1;
      if (n1 > 0) {
        i1++;
        win_data[i1][0] = (pre_index + 1) * win_size;
        win_data[i1][1] = win_data[i1][0] + n1 * win_size - 1;
        win_data[i1][2] = 0;
        win_data[i1][3] = n1;
      }

      i1++;
      win_data[i1][0] = j * win_size;
      win_data[i1][1] = win_data[i1][0] + win_size - 1;
      win_data[i1][2] = Weight[j] * 1.0 / win_size;
      win_data[i1][3] = 1;

      pre_index = j;
    }
  }
  num_win = i1;

}

void bcp_algo::buildChrMap(Reads& reads) {
  //treads should have been normalized against creads in the main app.cpp
  size_t chrNum = 0;
  foreach(std::string chr, reads.pos_reads.chrs()) {

    chr2Num[chr] = chrNum;
    num2Chr[chrNum] = chr;

    chrNum++;
  }
}

double bcp_algo::prob_pois(int l1, double lambda) {
  int i2;
  double prob, logprob, sum_temp;
  prob = exp(-1 * lambda);
  for (i2 = 1; i2 <= l1; i2++) {
    logprob = 0.0;
    sum_temp = st[i2];
    logprob = (-lambda) + i2 * log(lambda) - sum_temp;
    prob += exp(logprob);
  }
  return prob;
}

void bcp_algo::seg(Matrix<double> input, double thre) {
  int i, j, k, j1;
  double temp_sum;
  j = 1;
  i = 1;
  num_seg = 0;

  while (1) {
    temp_sum = 0.0;
    if (j > num_win)
      break;
    if (input[j][4] < thre) {
      j1 = j;
      while ((j1 <= num_win) && (input[j1][4] < thre))
        j1++;
      ss[i][0] = input[j][0];
      ss[i][1] = input[j1 - 1][1];
      ss[i][2] = ss[i][1] - ss[i][0] + 1.0;
      for (k = j; k < j1; k++)
        temp_sum += (input[k][1] - input[k][0] + 1.0) * input[k][4];
      ss[i][3] = temp_sum / ss[i][2];
      ss[i][4] = j;
      ss[i][5] = j1 - 1;
      i++;
      j = j1;
    } else {
      num_seg++;
      j1 = j;
      while ((j1 <= num_win) && (input[j1][4] >= thre))
        j1++;
      ss[i][0] = input[j][0];
      ss[i][1] = input[j1 - 1][1];
      ss[i][2] = ss[i][1] - ss[i][0] + 1.0;
      for (k = j; k < j1; k++)
        temp_sum += (input[k][1] - input[k][0] + 1.0) * input[k][4];
      ss[i][3] = temp_sum / ss[i][2];
      ss[i][4] = j;
      ss[i][5] = j1 - 1;
      i++;
      j = j1;

    }
  }
  num_allseg = i - 1;
}

int bcp_algo::frag_count1(int rec) {
  int start_ind, end_ind, comp_ind;
  int ind1, ind2;
  int count;

  start_ind = 1;
  end_ind = num_data_frag - 1;
  while (start_ind < end_ind - 1) {
    comp_ind = start_ind + (end_ind - start_ind) / 2;
    if (data_frag[comp_ind] == rseg[rec][1]) {
      start_ind = comp_ind;
      break;
    } else {
      if (data_frag[comp_ind] < rseg[rec][1])
        start_ind = comp_ind;
      else
        end_ind = comp_ind;
    }

  }
  ind1 = start_ind;

  start_ind = 1;
  end_ind = num_data_frag - 1;
  while (start_ind < end_ind - 1) {
    comp_ind = start_ind + (end_ind - start_ind) / 2;
    if ((data_frag[comp_ind] + frag_size - 1) == rseg[rec][0]) {
      end_ind = comp_ind;
      break;
    } else {
      if ((data_frag[comp_ind] + frag_size - 1) > rseg[rec][0])
        end_ind = comp_ind;
      else
        start_ind = comp_ind;
    }
  }
  ind2 = end_ind;

  count = ind1 - ind2 + 1;
  return count;
}

int bcp_algo::frag_count2(int rec) {
  int start_ind, end_ind, comp_ind;
  int ind1, ind2;
  int count;

  start_ind = 1;
  end_ind = num_input_frag - 1;
  while (start_ind < end_ind - 1) {
    comp_ind = start_ind + (end_ind - start_ind) / 2;
    if (input_frag[comp_ind] == rseg[rec][1]) {
      start_ind = comp_ind;
      break;
    } else {
      if (input_frag[comp_ind] < rseg[rec][1])
        start_ind = comp_ind;
      else
        end_ind = comp_ind;
    }

  }
  ind1 = start_ind;

  start_ind = 1;
  end_ind = num_input_frag - 1;
  while (start_ind < end_ind - 1) {
    comp_ind = start_ind + (end_ind - start_ind) / 2;
    if ((input_frag[comp_ind] + frag_size - 1) == rseg[rec][0]) {
      end_ind = comp_ind;
      break;
    } else {
      if ((input_frag[comp_ind] + frag_size - 1) > rseg[rec][0])
        end_ind = comp_ind;
      else
        start_ind = comp_ind;
    }
  }
  ind2 = end_ind;

  count = ind1 - ind2 + 1;
  return count;
}

bcp_algo::bcp_algo() :
  Weight(3000000, 0),frag_size(200),win_size(200),p_value(0.001),L1(3000000),L2(1500000),N(56) {
  }

bcp_algo::~bcp_algo() {
}

void bcp_algo::detectSummits(Reads& treatment_reads, Reads& control_reads,
    cmd_option_parser& option) {
  cmain(treatment_reads, control_reads, option);
}

