/*
 * @file conjugate_gradient_solver.h
 * @brief implementation of conjugate gradient solver for semi-definite sparse linear system Ax = b
 * @reference <An Introduction to the Conjugate Gradient Method without the Agonizing Pain>
 * @author Fei Zhu
 *
 * This file is part of Physika, a versatile physics simulation library.
 * Copyright (C) 2013- Physika Group.
 *
 * This Source Code Form is subject to the terms of the GNU General Public License v2.0.
 * If a copy of the GPL was not distributed with this file, you can obtain one at:
 * http://www.gnu.org/licenses/gpl-2.0.html
 *
 */

#ifndef PHYSIKA_NUMERICS_LINEAR_SYSTEM_SOLVERS_CONJUGATE_GRADIENT_SOLVER_H_
#define PHYSIKA_NUMERICS_LINEAR_SYSTEM_SOLVERS_CONJUGATE_GRADIENT_SOLVER_H_

#include "Physika_Numerics/Linear_System_Solvers/iterative_solver.h"

namespace Physika{

template <typename Scalar>
class ConjugateGradientSolver: public IterativeSolver<Scalar>
{
public:
    ConjugateGradientSolver();
    ~ConjugateGradientSolver();
    virtual bool solve(const LinearSystem<Scalar> &system, const GeneralizedVector<Scalar> &b, GeneralizedVector<Scalar> &x);
protected:
};

}  //end of namespace Physika

#endif //PHYSIKA_NUMERICS_LINEAR_SYSTEM_SOLVERS_CONJUGATE_GRADIENT_SOLVER_H_