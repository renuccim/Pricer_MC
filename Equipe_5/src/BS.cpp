#include "BS.h"
#include <iostream>

using namespace std;

BS::BS(Parser *P)
{
	P->extract("option size",this->size_);
	P->extract("interest rate",this->r_);
	P->extract("correlation",this->rho_);
	P->extract("volatility", this->sigma_, size_);	
	P->extract("spot", this->spot_, size_);
	P->extract("dividend rate", this->trend_, size_);
	// The cholesky factorization
	//L = pnl_mat_create_from_scalar(this->size_,this->size_,1);
	L = pnl_mat_create_from_scalar(size_, size_, rho_);

	for (int i = 0; i < size_; i++)
	{
		pnl_mat_set(L, i, i, 1);
	}
	pnl_mat_chol(L);
	
}

BS::BS(int size,double r,double rho,PnlVect *sigma,PnlVect *spot){
	
	size_ = size;
	r_ = r;
	rho_ = rho;
	sigma_ = sigma;
	spot_ = spot;

	L = pnl_mat_create_from_scalar(size_, size_, rho_);

	for (int i = 0; i < size_; i++)
	{
		pnl_mat_set(L, i, i, 1);
	}

}

BS::~BS()
{
#ifdef _DEBUG
	cout << "~BS() : Ready to call pnl_vect_free on sigma, spot and trend ->->-> " << endl;
#endif
	pnl_vect_free(&this->sigma_);
	pnl_vect_free(&this->spot_);
	pnl_vect_free(&this->trend_);	
	pnl_mat_free(&L);
#ifdef _DEBUG
	cout << "~BS() : Successfull call of pnl_rng_free" << endl;
#endif
}


void BS::asset(PnlMat *path, double T, int N, PnlRng *rng){

	//cout << "in BS" << endl;
	
	//nb rows and columns unreachable

	pnl_mat_resize(path, N+1, size_);
	pnl_mat_set_row(path, spot_, 0);

	PnlVect *W = pnl_vect_create(size_);

	double Sti_1;
	double e;
	double sigma_d;
	double t_  = T/double(N);

	double LdW;

	PnlVect *Ld = pnl_vect_create(size_);
	//cout << "before BS" << endl;
	for (int i = 1; i < N+1; i++)
	{

		//cout << "A l instant i = " << i <<"  ";
		//cout << "before rng BS" << endl;
		pnl_vect_rng_normal(W,size_,rng);
		//pnl_vect_print(W);
		//cout << "after rng BS" << endl;
		
		for (int d = 0; d < size_; d++)
		{

			//cout << " d = " << d ;
			Sti_1 = pnl_mat_get(path, i-1, d);
			//cout << "   Sti_1_  : " << Sti_1 ;
			//cout << "get vect BS" << endl;
			sigma_d = pnl_vect_get(sigma_,d);
			//cout << "   sigma_ : " << sigma_d;
			pnl_mat_get_row(Ld,L,d);
			//cout << "scalar prod BS" << endl;
			LdW = pnl_vect_scalar_prod(Ld,W);
			//cout << "   LdW : " << LdW;		
			e = exp((r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW);
			//cout << "   Valeur " << exp(r_ - 0.5 * pow(sigma_d,2.0) ) * t_ ;
			//cout << "      ";
			pnl_mat_set(path, i, d, Sti_1*e);
		}
		//cout << "" << endl;
	}
	//cout << "in out BS" << endl;
	//pnl_vect_free(&V);
	pnl_vect_free(&Ld);
	pnl_vect_free(&W);


}

void BS::asset(PnlMat *path, double t, int N, double T,
          PnlRng *rng, const PnlMat *past){

	pnl_mat_resize(path, size_, N+1);
	pnl_mat_set_subblock(path, past, int(N*t), size_);


	//PnlMat *L = pnl_mat_create_from_scalar(size_, size_, rho_);
	PnlVect *W = pnl_vect_create(size_);

	// for (int i = 0; i < size_; i++)
	// {
	// 	pnl_mat_set(L, i, i, 1);
	// }

	pnl_mat_chol(L);

	double Sti_1;
	double e;
	double sigma_d;
	double LdW;
	PnlVect *Ld;

	int i = int(N*t)+1;
	double t_  = T * (double(i)/double(N)) - t;

	pnl_vect_rng_normal(W,size_,rng);

	for (int d = 0; d < size_; d++)
	{
		Sti_1 = pnl_mat_get(past, int(N*t), d);
		sigma_d = pnl_vect_get(sigma_,d);
		pnl_mat_get_row(Ld,L,d);
		LdW = pnl_vect_scalar_prod(Ld,W);
		e = exp(r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW;
		pnl_mat_set(path, i, d, Sti_1*e);
	}

	t_  = T/double(N);

	for (int i = int(N*t) +2 ; i < N+1; i++)
	{
		pnl_vect_rng_normal(W,size_,rng);

		for (int d = 0; d < size_; d++)
		{
			Sti_1 = pnl_mat_get(path, i-1, d);
			sigma_d = pnl_vect_get(sigma_,d);
			pnl_mat_get_row(Ld,L,d);
			LdW = pnl_vect_scalar_prod(Ld,W);
			e = exp(r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW;
			pnl_mat_set(path, i, d, Sti_1*e);
		}
	}

	//pnl_mat_free(L);
	//pnl_vect_free(&V);
	pnl_vect_free(&Ld);
	pnl_vect_free(&W);
	
}


void BS::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *ic){

	// int size = mod_.size_;
	// int N = opt_.TimeSteps_;

	// PnlMat *path = pnl_mat_create(size,N);

	// for (int i = 0; i < samples_; i++)
	// {
	// 	mod_.asset(path, t, N, T, rng, past);
	// 	for (int d = 0; d < size; d++)
	// 	{
	// 		for (int k = int(N*t)+1; k < ; k++)
	// 		{
	// 			pnl_mat_set(path, k, d, (1+h_) * pnl_mat_get(path, d, k));
	// 		}

	// 		pnl_vect_set(delta, d, opt_.payoff(path) + pnl_vect_get(delta,d));
			
	// 		for (int k = int(N*t)+1; k < ; k++)
	// 		{
	// 			pnl_mat_set(path, k, d, (1-h_)/(1+h_) * pnl_mat_get(path, d, k));
	// 		}
			
	// 		pnl_vect_set(delta, d, pnl_vect_get(delta,d) - opt_.payoff(path) );
	// 	}

	// }

	// for (int d = 0; d < size; d++)
	// {
	// 	pnl_vect_set(delta, d, exp((-mod_.r_)*(opt_.T_-t))/(samples_* 2 * h_ * pnl_mat_get(int(N*t),d)));

	// }

	// pnl_mat_free(path);
}

void BS::shift_asset(PnlMat *shift_path, const PnlMat *path,
                   int d, double h, double t, double timestep){

}