#include <iostream>
#include "parser.h"
#include "BS.h"


using namespace std;

int main(int argc, char **argv)
{
  char *infile = argv[1];
  Parser *P = new Parser(infile);

  BS *mod_ = new BS(P);

  cout << "option size " << mod_->size_ << endl;
  cout << "interest rate " << mod_->r_ << endl;
  cout << "correlation " << mod_->rho_ << endl;
  cout << "spot "; pnl_vect_print_asrow(mod_->spot_);
  cout << "volatility "; pnl_vect_print_asrow(mod_->sigma_);
  cout << "L matrix " << endl; pnl_mat_print(mod_->L);

  int N = 10;

  PnlMat *path = pnl_mat_create(N,mod_->size_);
  PnlRng *rng;// = pnl_rng_new();
  double T;

  P->extract("maturity", T);

  mod_->asset(path, T, N, rng);

  pnl_mat_print(path);

  pnl_mat_free(&path);
  pnl_rng_free(&rng);
  delete P;
  delete mod_;
  exit(0);
}
