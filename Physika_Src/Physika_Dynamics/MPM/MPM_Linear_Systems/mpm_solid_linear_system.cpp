/*
 * @file mpm_solid_linear_system.cpp
 * @brief linear system for implicit integration of MPMSolid && CPDIMPMSolid driver
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

#include <typeinfo>
#include <vector>
#include "Physika_Core/Matrices/matrix_2x2.h"
#include "Physika_Core/Matrices/matrix_3x3.h"
#include "Physika_Core/Utilities/physika_exception.h"
#include "Physika_Dynamics/Particles/solid_particle.h"
#include "Physika_Dynamics/MPM/mpm_solid.h"
#include "Physika_Dynamics/MPM/MPM_Linear_Systems/mpm_uniform_grid_generalized_vector.h"
#include "Physika_Dynamics/MPM/MPM_Linear_Systems/mpm_solid_linear_system.h"

namespace Physika{

template <typename Scalar, int Dim>
MPMSolidLinearSystem<Scalar,Dim>::MPMSolidLinearSystem(MPMSolid<Scalar,Dim> *mpm_solid_driver)
    :mpm_solid_driver_(mpm_solid_driver), active_obj_idx_(-1)
{
}

template <typename Scalar, int Dim>
MPMSolidLinearSystem<Scalar,Dim>::~MPMSolidLinearSystem()
{

}

template <typename Scalar, int Dim>
void MPMSolidLinearSystem<Scalar,Dim>::multiply(const GeneralizedVector<Scalar> &x,
                                                GeneralizedVector<Scalar> &result) const
{
    try
    {
        const MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> > &xx = dynamic_cast<const MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> >&>(x);
        MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> > &rr = dynamic_cast<MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> >&>(result);
        energyHessianMultiply(xx, rr);
        Scalar dt_square = mpm_solid_driver_->computeTimeStep();
        dt_square *= dt_square;
        std::vector<Vector<unsigned int, Dim> > active_grid_nodes;
        if (active_obj_idx_ == -1) //all objects solved together
        {
            mpm_solid_driver_->activeGridNodes(active_grid_nodes);
            for (unsigned int i = 0; i < active_grid_nodes.size(); ++i)
            {
                Vector<unsigned int, Dim> &node_idx = active_grid_nodes[i];
                rr[node_idx] = xx[node_idx] + dt_square*rr[node_idx] / mpm_solid_driver_->gridMass(node_idx);
            }
        }
        else  //solve for active object
        {
            mpm_solid_driver_->activeGridNodes(active_obj_idx_, active_grid_nodes);
            for (unsigned int i = 0; i < active_grid_nodes.size(); ++i)
            {
                Vector<unsigned int, Dim> &node_idx = active_grid_nodes[i];
                rr[node_idx] = xx[node_idx] + dt_square*rr[node_idx] / mpm_solid_driver_->gridMass(active_obj_idx_, node_idx);
            }
        }
    }
    catch (std::bad_cast &e)
    {
        throw PhysikaException("Incorrect argument!");
    }
}

template <typename Scalar, int Dim>
void MPMSolidLinearSystem<Scalar,Dim>::preconditionerMultiply(const GeneralizedVector<Scalar> &x,
                                                              GeneralizedVector<Scalar> &result) const
{
    try{
        const MPMUniformGridGeneralizedVector<Vector<Scalar,Dim> > &mpm_x = dynamic_cast<const MPMUniformGridGeneralizedVector<Vector<Scalar,Dim> >&>(x);
        MPMUniformGridGeneralizedVector<Vector<Scalar,Dim> > &mpm_result = dynamic_cast<MPMUniformGridGeneralizedVector<Vector<Scalar,Dim> >&>(result);
        jacobiPreconditionerMultiply(mpm_x,mpm_result);
    }
    catch(std::bad_cast &e)
    {
        throw PhysikaException("Incorrect argument type!");
    }
}

template <typename Scalar, int Dim>
void MPMSolidLinearSystem<Scalar,Dim>::setActiveObject(int obj_idx)
{
    active_obj_idx_ = obj_idx;
    //all negative values are set to -1
    active_obj_idx_ = active_obj_idx_ < 0 ? -1 : active_obj_idx_;
}

template <typename Scalar, int Dim>
void MPMSolidLinearSystem<Scalar, Dim>::energyHessianMultiply(const MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> > &x_diff,
                                                                  MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> > &result) const
{
    std::vector<Vector<unsigned int, Dim> > active_grid_nodes;
    std::vector<SquareMatrix<Scalar, Dim> > A_p(mpm_solid_driver_->totalParticleNum(), SquareMatrix<Scalar, Dim>(0));
    unsigned global_particle_idx = 0;
    std::vector<Vector<unsigned int, Dim> > nodes_in_range;
    for (unsigned int obj_idx = 0; obj_idx < mpm_solid_driver_->objectNum(); ++obj_idx)
    {
        for (unsigned int particle_idx = 0; particle_idx < mpm_solid_driver_->particleNumOfObject(obj_idx); ++particle_idx, ++global_particle_idx)
        {
            mpm_solid_driver_->gridNodesInRange(obj_idx, particle_idx,nodes_in_range);
        }
    }
    if (active_obj_idx_ == -1) //all objects solved together
    {
        mpm_solid_driver_->activeGridNodes(active_grid_nodes);
    }
    else  //solve for one active object
    {
        mpm_solid_driver_->activeGridNodes(active_obj_idx_, active_grid_nodes);
    }
}

template <typename Scalar, int Dim>
void MPMSolidLinearSystem<Scalar,Dim>::jacobiPreconditionerMultiply(const MPMUniformGridGeneralizedVector<Vector<Scalar,Dim> > &x,
                                                                    MPMUniformGridGeneralizedVector<Vector<Scalar, Dim> > &result) const
{
    if (active_obj_idx_ == -1) //all objects solved together
    {
        //TO DO
    }
    else  //solve for one active object
    {
        //TO DO
    }
}

//explicit instantiations
template class MPMSolidLinearSystem<float,2>;
template class MPMSolidLinearSystem<float,3>;
template class MPMSolidLinearSystem<double,2>;
template class MPMSolidLinearSystem<double,3>;

}  //end of namespace Physika