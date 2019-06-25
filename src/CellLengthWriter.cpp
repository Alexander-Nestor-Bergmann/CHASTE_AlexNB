/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellLengthWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
// #include <math.h>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellLengthWriter<ELEMENT_DIM, SPACE_DIM>::CellLengthWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllength.dat")
{
    this->mVtkCellDataName = "Cell Length";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CellLengthWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutputAsVector(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // We first need to test that the cell population is vertex-based
    VertexBasedCellPopulation<SPACE_DIM>* p_vbcp = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(pCellPopulation);

    if (p_vbcp == nullptr)
    {
        EXCEPTION("Cell length is only associated with vertex-based cell populations");
    }

    // Get the short axis of the element
    AbstractMesh<SPACE_DIM, SPACE_DIM>& r_mesh = p_vbcp->rGetMesh();
    // Pointer to mesh
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* p_mesh = static_cast<MutableVertexMesh<SPACE_DIM, SPACE_DIM>*>(&r_mesh);
    // Toroidal2dVertexMeshWithMutableSize* p_mesh = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);

    // Get the element corresponding to the cell
    unsigned index = p_vbcp->GetLocationIndexUsingCell(pCell);



    c_vector<double, SPACE_DIM> short_axis = zero_vector<double>(SPACE_DIM);

    // Calculate the moments of the element about its centroid (recall that I_xx and I_yy must be non-negative)
    c_vector<double, 3> moments = p_mesh->CalculateMomentsOfElement(index);

    // Normalise the moments vector to remove problem of a very small discriminant (see #2874)
    moments /= norm_2(moments);

    // If the principal moments are equal...
    double discriminant = (moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2);
    if (fabs(discriminant) < DBL_EPSILON)
    {
        // ...then every axis through the centroid is a principal axis, so return a random unit vector
        short_axis(0) = 1;
        short_axis(1) = 0;
    }
    else
    {
        // If the product of inertia is zero, then the coordinate axes are the principal axes
        if (fabs(moments(2)) < DBL_EPSILON)
        {
            if (moments(0) < moments(1))
            {
                short_axis(0) = 0.0;
                short_axis(1) = 1.0;
            }
            else
            {
                short_axis(0) = 1.0;
                short_axis(1) = 0.0;
            }
        }
        else
        {
            // Otherwise we find the eigenvector of the inertia matrix corresponding to the largest eigenvalue
            double lambda = 0.5 * (moments(0) + moments(1) + sqrt(discriminant));

            short_axis(0) = 1.0;
            short_axis(1) = (moments(0) - lambda) / moments(2);

            // Normalise the short axis before returning it
            // short_axis /= norm_2(short_axis);
        }
    }


    // // Get the short axis
    // c_vector<double, SPACE_DIM> short_axis = p_mesh->GetShortAxisOfElement(elem_index);

    // Dot product with AP and Dv directions
    c_vector<double, SPACE_DIM> AP_DV_directions;
    AP_DV_directions[0] = -short_axis[0]; // AP
    AP_DV_directions[1] = short_axis[1]; // DV

    // unsigned cell_length = p_vbcp->GetRosetteRankOfCell(pCell);

    return AP_DV_directions;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellLengthWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);

    // c_vector<double, SPACE_DIM> cell_lengths = this->GetCellDataForVtkOutput(pCell, pCellPopulation);
    c_vector<double, SPACE_DIM> cell_lengths = this->GetCellDataForVtkOutputAsVector(pCell, pCellPopulation);


    *this->mpOutStream << location_index << " " << cell_id << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << centre_location[i] << " ";
    }

    *this->mpOutStream << cell_lengths[0] << " ";
    *this->mpOutStream << cell_lengths[1] << " ";
    // *this->mpOutStream << cell_lengths << " ";
}

// Explicit instantiation
template class CellLengthWriter<1,1>;
template class CellLengthWriter<1,2>;
template class CellLengthWriter<2,2>;
template class CellLengthWriter<1,3>;
template class CellLengthWriter<2,3>;
template class CellLengthWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellLengthWriter)
