#ifndef _MC_H
#define _MC_H

#include "option.h"
#include "asian.h"
#include "basket.h"
#include "barrier_l.h"
#include "barrier_u.h"
#include "barrier.h"
#include "performance.h"
#include "bs.h"
#include "parser.h"
#include <cmath>
#include <cassert>
#include <ctime>
#include <iomanip> 
#include "pnl/pnl_random.h"
#include <cassert>

class MonteCarlo
{
public:
  BS *mod_; /*! pointeur vers le mod�le */
  Option *opt_; /*! pointeur sur l'option */
  PnlRng *rng; /*! pointeur sur le g�n�rateur */
  double h_; /*! pas de diff�rence finie */
  int samples_; /*! nombre de tirages Monte Carlo */

	MonteCarlo(Parser *P);
	~MonteCarlo();

  /**
   * Calcule le prix de l'option � la date 0
   *
   * @param[out] prix valeur de l'estimateur Monte Carlo
   * @param[out] ic largeur de l'intervalle de confiance
   */
  void price(double &prix, double &ic);

  /**
   * Calcule le prix de l'option � la date t
   *
   * @param[in]  past contient la trajectoire du sous-jacent
   * jusqu'� l'instant t
   * @param[in] t date � laquelle le calcul est fait
   * @param[out] prix contient le prix
   * @param[out] ic contient la largeur de l'intervalle
   * de confiance sur le calcul du prix
   */
  void price(const PnlMat *past, double t, double &prix, double &ic);


  /**
   * Calcule le delta de l'option � la date t
   *
   * @param[in] past contient la trajectoire du sous-jacent
   * jusqu'� l'instant t
   * @param[in] t date � laquelle le calcul est fait
   * @param[out] delta contient le vecteur de delta
   */
  void delta_(const PnlMat *past, double t, PnlVect *delta);

	/**
   * Calcule le delta de l'option � la date t
   *
   * @param[in] past contient la trajectoire du sous-jacent
   * jusqu'� l'instant t
   * @param[in] t date � laquelle le calcul est fait
   * @param[out] delta contient le vecteur de delta
   */
  void delta(const PnlMat *past, double t, PnlVect *delta);

	/**
   * Construit le portefeuille de couverture et calcule le P&L
   *
	 * @param[in]  H nombre de date de constatation
	 * @param[in]  marketPath matrice de taille d x (H+1) qui contient
	 * une simulation du march�
   * @param[out] V vecteur des valeurs de portefeuille de couverture
	 * de dimension H+1
   * @param[out] PL Profit and Loss
   */
  void hedge(PnlVect *V, double &PL, int H, const PnlMat *marketPath);

	/**
   * Construit le portefeuille de couverture et calcule le P&L
   *
	 * @param[in]  H nombre de date de constatation
   * @param[out] V vecteur des valeurs de portefeuille de couverture
	 * de dimension H+1
   * @param[out] PL Profit and Loss
   */
  void hedge(PnlVect *V, double &PL, int H);
};

#endif /* _MC_H */

