#include <iostream>
#include <time->h>
#include "BS.h"

using namespace std;

BS::BS(Parser myParser)
{
	P->extract("option size",this->size_);
	P->extract("interest rate",this->r_);
	P->extract("correlation",this->rho_);
	P->extract("sigma",this->sigma_, this->size_);	
	P->extract("spot",this->spot_, this->size_);
	P->extract("trend",this->trend_, this->size_);
	// The cholesky factorization
	L = pnl_mat_create_from_scalar(this->size_,this->size_,1);
	for(int i=1; i <= this->size_; i++)
	{
		for(int j=1; j <= this->size_; j++
		{
			if (i!=j)
				MLET(Correlation,i,j) = this->rho_;
		}
	}
	pnl_mat_chol(L);
	
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

void BS::asset(PnlMat *path, double T, int N, PnlRng *rng)
{
	double step = T/N;
	double prodScal = 0;
	double sigma_d = 0;

	// The Gaussian vector
	PnlVect *G = pnl_vect_create_from_zero(this->size_);
	PnlVect *Ld = pnl_vect_create_from_zero(this->size_);	
	// First col of path contains spot
	for(int d=1; d <= this->size_; d++)
		MLET(path,d,1) = GET(this->spot_,d);

	// Generation from 2 to N+1(th) column
	for(int ti=2; ti <= N+1; ti++)
	{
		// Gaussian Dimension size_ generation
		pnl_vect_rng_normal(G,this->size_,rng);
		for(int d=1; d <= this->size_; d++)
		{
			pnl_mat_get_col(Ld,L,d);
			prodScal = pnl_vect_scalar_prod(Ld,G)
			sigma_d = GET(this->spot_,d);
			MLET(path,d,ti) = MGET(path,d,ti-1)*exp( (this->r_-pow(sigma_d,2)/2)*step + sigma_d*sqrt(step)*prodScal );
		}
	}
	
	// Memory free
	pnl_vect_free(&G);
	pnl_vect_free(&Ld);
}
