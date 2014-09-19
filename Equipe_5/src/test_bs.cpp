#include <iostream>
#include "parser.h"
#include "bs.h"
#include "option.h"
#include "asian.h"
#include "basket.h"
#include "barrier_l.h"
#include "barrier_u.h"
#include "barrier.h"
#include "performance.h"


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
  PnlRng *rng = pnl_rng_create(PNL_RNG_MERSENNE);
  pnl_rng_sseed(rng,0);

  double T;

  P->extract("maturity", T);

  PnlMat *past = pnl_mat_create(1,mod_->size_);
  pnl_mat_set_row(past,mod_->spot_,0);
  cout << "" << endl;
  pnl_mat_print(past);
  mod_->asset(path,0.0, N, T, rng, past);
  pnl_mat_print(path);


  pnl_mat_free(&past);
  pnl_mat_free(&path);
  pnl_rng_free(&rng);
  delete P;
  delete mod_;
exit(0);
}
