#ifndef _BARRIER_U_H
#define _BARRIER_U_H

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include <cmath>
#include "parser.h"
#include "option.h"

/// \brief Classe Option BARRIER_U
class Barrier_u : public Option
{
public:
  double K_; /// Strike
  PnlVect *payoffCoefficients_; /// Coefficients intervenant dans le calcul du payoff
	PnlVect *upperBarrier_; /// Upper barrier

	Barrier_u(Parser *P);
	~Barrier_u();

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


#endif /* _BARRIER_U_H */
