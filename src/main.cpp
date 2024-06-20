#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <random>
#include "doe_optimizer.h"

#include <Rcpp.h>
#include <R.h>

using namespace Rcpp;
using namespace std;


int criteria_selector(string crit)
{
  int critopt;
  critopt = 1;
  if(crit == "CD2") critopt=1;
  else if(crit == "MD2") critopt=2;
  else if(crit == "WD2") critopt=3;
  else if(crit == "maximin") critopt=4;
  else if(crit == "MC") critopt=5;
  else if(crit == "A2") critopt=6;
  return critopt;
}

NumericMatrix Generate_init_matrix(string init_method, int nsamp, int nv, int nlevel, NumericMatrix initX, int rand_seed)
{
  int i,j;
  vector<double> col;
  NumericMatrix return_matrix = NumericMatrix(nsamp,nv);

  if (init_method == "input" && initX.nrow()>1)
  {
    for(i=0;i<initX.ncol();i++) {
      for(j=0;j<initX.nrow();j++) {
        return_matrix(j,i) = initX(j,i);
      }
    }
  }
  else if (init_method == "rand")
  {
    for(i=1;i<=nsamp;i++) col.push_back((i%nlevel)+1);
    for(i=0;i<nv;i++)
    {
      std::mt19937 rng(rand_seed);
      std::shuffle(col.begin(), col.end(), rng);
      for(j=0;j<nsamp;j++)  return_matrix(j,i) = col[j];
    }
  }
  return return_matrix;
}

NumericMatrix Generate_Aug_matrix(string init_method, NumericMatrix xp, int nnew, int nv,
                                  int nlevel, NumericMatrix initX)
{
    int i,j,k, np, nsamp, max_k=0, fill_size, all_fill_size;
    np = xp.nrow();
    nsamp = np + nnew;
    vector<double> temp(np,0);

    vector<vector<int> > freq_table;
    NumericMatrix return_matrix  = NumericMatrix(nnew,nv);
    freq_table.assign(nlevel, vector<int>(nv, 0));

    if (init_method == "input" && initX.nrow()>1)
    {
      for(i=0;i<initX.ncol();i++) {
        for(j=0;j<initX.nrow();j++) {
          return_matrix(j,i) = initX(j,i);
        }
      }
    }
    else if (init_method == "rand")
    {
      for (i=0;i<nv;i++)
      {
        for (j=0;j<np;j++) temp[j] = xp(j,i);
        for (j=0;j<nlevel;j++)
        {
          freq_table[j][i] = (int) count(temp.begin(), temp.end(), j+1);
          if (max_k<freq_table[j][i]) max_k = freq_table[j][i];
        }
      }
      max_k = max(max_k, nsamp/nlevel);
      for (i=0;i<nv;i++)
      {
        all_fill_size = 0;
        for (j=0;j<nlevel;j++)
        {
          fill_size = max_k-freq_table[j][i];
          for (k=all_fill_size;k<all_fill_size+fill_size;k++)
          {
            return_matrix(k,i) = j+1;
          }
          all_fill_size += fill_size;
        }
      }
    }
    return return_matrix;
}

// [[Rcpp::export]]
double CritEval(NumericMatrix x0, int nlevel, StringVector crit)
{
  Criteria *c;
  double criteria = 1.0;
  int i, j, nv = x0.ncol(), nsamp= x0.nrow();
  int critopt = criteria_selector(as<string> (crit));
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));

  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) x[i][j] = x0(i,j);
  switch(critopt)
  {
  case 1:
    c = new CD2(x, nsamp, nv, nlevel);
    break;
  case 2:
    c = new MD2(x, nsamp, nv, nlevel);
    break;
  case 3:
    c = new WD2(x, nsamp, nv, nlevel);
    break;
  case 4:
    c = new Maximin(x, nsamp, nv, nlevel);
    break;
  case 5:
    c = new MC(x, nsamp, nv, nlevel);
    break;
  case 6:
    c = new A2(x, nsamp, nv, nlevel);
    break;
  default:
    c = new CD2(x, nsamp, nv, nlevel);
  break;
  }
  criteria = c->get_criteria();
  x.clear();
  x.shrink_to_fit();
  delete c;
  return(criteria);
}

// [[Rcpp::export]]
List SATA_UD(int nsamp, int nv, int nlevel, StringVector init_method, NumericMatrix initX,
             StringVector crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed)
{
  int i,j;
  List lst;
  clock_t start_time;
  vector<double> critobj_vector;
  vector<int> optimize_columns(nv,1);
  vector<vector<double> > final_design;
  double critobj, critobj0, search_time;

  int critopt = criteria_selector(as<string> (crit));
  NumericMatrix Init_matrix,return_matrix = NumericMatrix(nsamp,nv);
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));

  srand(rand_seed);
  start_time = clock();
  Init_matrix = Generate_init_matrix(as<string>(init_method),nsamp,nv,nlevel,initX,rand_seed);
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) x[i][j] = Init_matrix(i,j);
  Optimizer opt(x, nsamp, 0, nv, nlevel, optimize_columns, critopt, maxiter, hits_ratio, levelpermt);
  critobj_vector = opt.SATA_Optimize();
  final_design = opt.get_design();
  critobj0 = critobj_vector.front();
  critobj = critobj_vector.back();
  for(i=0;i<nsamp; i++) for(j=0;j<nv;j++) return_matrix(i,j) = final_design[i][j];

  search_time = (double)(clock()-start_time)/CLOCKS_PER_SEC;

  lst["Init_Design"] = Init_matrix;
  lst["Final_Design"] = return_matrix;
  lst["Init_Obj"] = critobj0;
  lst["Final_Obj"] = critobj;
  lst["Time_Second"]= search_time;
  lst["Criterion_history"] = wrap(critobj_vector);
  return lst;
}

// [[Rcpp::export]]
List SATA_AUD(NumericMatrix xp,int nnew, int nv, int nlevel, StringVector init_method, NumericMatrix initX,
              StringVector crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed)
{
  int i,j;
  List lst;
  int np=xp.nrow();
  int nsamp = xp.nrow()+nnew;
  clock_t start_time;
  vector<double> critobj_vector;
  vector<int> optimize_columns(nv,1);
  vector<vector<double> > final_design;
  double critobj, critobj0, search_time;
  int critopt = criteria_selector(as<string> (crit));
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));
  NumericMatrix InputX(nnew,nv),Init_matrix(nsamp,nv),return_matrix = NumericMatrix(nsamp,nv);

  srand(rand_seed);
  start_time = clock();
  InputX = Generate_Aug_matrix(as<string>(init_method),xp,nnew,nv,nlevel,initX);
  for(j=0;j<nv;j++)
  {
    for(i=0;i<np;i++) {x[i][j] =Init_matrix(i,j)= xp(i,j);}
    for(i=0;i<nnew;i++) {x[i+np][j] = Init_matrix(i+np,j)=InputX(i,j);}
  }
  Optimizer opt(x, nnew, np, nv, nlevel, optimize_columns, critopt, maxiter, hits_ratio, levelpermt);
  critobj_vector = opt.SATA_Optimize();
  final_design = opt.get_design();
  critobj0 = critobj_vector[0];
  critobj = critobj_vector.back();
  for(i=0;i<nsamp; i++) for(j=0;j<nv;j++) return_matrix(i,j) = final_design[i][j];
  search_time = (double)(clock()-start_time)/CLOCKS_PER_SEC;

  lst["Init_Design"] = Init_matrix;
  lst["Final_Design"] = return_matrix;
  lst["Init_Obj"] = critobj0;
  lst["Final_Obj"] = critobj;
  lst["Time_Second"]= search_time;
  lst["Criterion_history"] = wrap(critobj_vector);
  return lst;
}

// [[Rcpp::export]]
List SATA_AUD_COL(NumericMatrix xp, int nvnew, int nlevel, StringVector init_method, NumericMatrix initX,
                  StringVector crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed)
{
  int i,j;
  List lst;
  int nsamp = xp.nrow();
  int nvp = xp.ncol();
  int nv = nvnew+nvp;
  clock_t start_time;
  vector<double> critobj_vector;
  vector<int> optimize_columns(nv,1);
  vector<vector<double> > final_design;
  double critobj, critobj0, search_time;
  int critopt = criteria_selector(as<string> (crit));
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));
  NumericMatrix InputX(nsamp,nvnew),Init_matrix(nsamp,nv),return_matrix = NumericMatrix(nsamp,nv);

  srand(rand_seed);
  start_time = clock();
  InputX = Generate_init_matrix(as<string>(init_method),nsamp,nvnew,nlevel,initX,rand_seed);
  for(j=0;j<nvp;j++) optimize_columns[j] = 0;
  for(i=0;i<nsamp;i++)
  {
    for(j=0;j<nvp;j++) x[i][j] =Init_matrix(i,j)= xp(i,j);
    for(j=0;j<nvnew;j++) x[i][j+nvp] =Init_matrix(i,j+nvp)= InputX(i,j);
  }

  Optimizer opt(x, nsamp, 0, nv, nlevel, optimize_columns, critopt, maxiter, hits_ratio, levelpermt);
  critobj_vector = opt.SATA_Optimize();
  final_design = opt.get_design();
  critobj0 = critobj_vector[0];
  critobj = critobj_vector.back();
  for(i=0;i<nsamp; i++) for(j=0;j<nv;j++) return_matrix(i,j) = final_design[i][j];
  search_time = (double)(clock()-start_time)/CLOCKS_PER_SEC;

  lst["Init_Design"] = Init_matrix;
  lst["Final_Design"] = return_matrix;
  lst["Init_Obj"] = critobj0;
  lst["Final_Obj"] = critobj;
  lst["Time_Second"]= search_time;
  lst["Criterion_history"] = wrap(critobj_vector);
  return lst;
}
