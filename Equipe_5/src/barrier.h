#ifndef _BARRIER_H
#define _BARRIER_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cmath>
#include "parser.h"
#include "option.h"
#include <cassert>

/// \brief Classe Option BARRIER
class Barrier : public Option
{
public:
  double K_; /// Strike
  PnlVect *payoffCoefficients_; /// Coefficients intervenant dans le calcul du payoff
	PnlVect *lowerBarrier_; /// Lower barrier
	PnlVect *upperBarrier_; /// Upper barrier

	Barrier(Parser *P);
	virtual ~Barrier();

  /**
   * Calcule la valeur du payoff sur la trajectoire
   *
   * @param[in] path est une matrice de taille d x (N+1)
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * @return phi(trajectoire)
   */
  double payoff(const PnlMat *path);
	
	void print();
};


#endif /* _BARRIER_H */
