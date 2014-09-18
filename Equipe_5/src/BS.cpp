#include "BS.h"

sing namespace std;


void BS::asset(PnlMat *path, double T, int N, PnlRng *rng){
	
	//nb rows and columns unreachable

	pnl_mat_resize(path, size_, N+1);
	pnl_mat_set_col(path, spot_, 0);

	PnlMat *L = pnl_mat_create_from_scalar(size_, N, rho_);
	PnlVect *W = pnl_vect_create(_size);
	//pnl_vect_rng_normal(W,size_,rng);

	for (int i = 0; i < size_; i++)
	{
		pnl_mat_set(L, i, i, 1);
		//pnl_vect_set(W,i,pnl_rng_normal(rng));
	}

	pnl_mat_chol(L);

	double Sti_1;
	double e;
	double sigma_d;
	double t_  = T/double(N);
	double LdW;

	for (int i = 1; i < N+1; i++)
	{
		pnl_vect_rng_normal(W,size_,rng);
		
		for (int d = 0; d < size_; d++)
		{
			Sti_1 = pnl_mat_get(path, i-1, d);
			sigma_d = pnl_vect_get(sigma_,d);
			LdW = pnl_vect_scalar_prod(pnl_vect_wrap_mat_row(L,d),W);
			e = exp(r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW;
			pnl_mat_set(path, i, d, Sti_1*e);
		}
	}

	pnl_mat_free(L);
	pnl_vect_free(W);


}

void BS::asset(PnlMat *path, double t, int N, double T,
          PnlRng *rng, const PnlMat *past){

	pnl_mat_resize(path, size_, N+1);
	pnl_mat_set_subblock(path, past, _size, int(N*t))


	PnlMat *L = pnl_mat_create_from_scalar(size_, N, rho_);
	PnlVect *W = pnl_vect_create(_size);

	for (int i = 0; i < size_; i++)
	{
		pnl_mat_set(L, i, i, 1);
	}

	pnl_mat_chol(L);

	double Sti_1;
	double e;
	double sigma_d;
	double LdW;

	int i = int(N*t)+1;
	double t_  = T * (double(i)/double(N)) - t;

	pnl_vect_rng_normal(W,size_,rng);

	for (int d = 0; d < size_; d++)
	{
		Sti_1 = pnl_mat_get(path, i-1, d);
		sigma_d = pnl_vect_get(sigma_,d);
		LdW = pnl_vect_scalar_prod(pnl_vect_wrap_mat_row(L,d),W);
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
			LdW = pnl_vect_scalar_prod(pnl_vect_wrap_mat_row(L,d),W);
			e = exp(r_ - 0.5 * pow(sigma_d,2.0) ) * t_ + sigma_d * sqrt(t_) * LdW;
			pnl_mat_set(path, i, d, Sti_1*e);
		}
	}

	pnl_mat_free(L);
	pnl_vect_free(W);
	
}

void BS::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *ic){

}

void BS::shift_asset(PnlMat *shift_path, const PnlMat *path,
                   int d, double h, double t, double timestep){

}