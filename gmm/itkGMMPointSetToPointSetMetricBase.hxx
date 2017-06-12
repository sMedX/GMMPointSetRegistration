/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkGMMPointSetToPointSetMetricBase_hxx
#define itkGMMPointSetToPointSetMetricBase_hxx

#include "itkGMMPointSetToPointSetMetricBase.h"

namespace itk
{
/** Constructor */
template< typename TFixedPointSet, typename TMovingPointSet >
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::GMMPointSetToPointSetMetricBase()
{
  m_FixedPointSet  = ITK_NULLPTR;    // has to be provided by the user.
  m_MovingPointSet = ITK_NULLPTR;    // has to be provided by the user.
  m_Transform      = ITK_NULLPTR;    // has to be provided by the user.

  m_Jacobian.set_size(MovingPointSetDimension, m_NumberOfParameters);
  m_JacobianCache.set_size(MovingPointSetDimension, MovingPointSetDimension);
}

/** Set the parameters that define a unique transform */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::SetTransformParameters(const ParametersType & parameters) const
{
  if ( !m_Transform )
    {
    itkExceptionMacro(<< "Transform has not been assigned");
    }
  m_Transform->SetParameters(parameters);
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::Initialize(void)
throw ( ExceptionObject )
{
  if (!m_Transform) {
    itkExceptionMacro(<< "Transform is not present");
  }

  m_NumberOfParameters = m_Transform->GetNumberOfParameters();

  if ( !m_MovingPointSet )
    {
    itkExceptionMacro(<< "MovingPointSet is not present");
    }

  if ( !m_FixedPointSet )
    {
    itkExceptionMacro(<< "FixedPointSet is not present");
    }

  // If the point set is provided by a source, update the source.
  if ( m_FixedPointSet->GetSource() )
    {
    m_FixedPointSet->GetSource()->Update();
    }

  m_FixedPointMatrix.set_size(m_FixedPointSet->GetNumberOfPoints(), FixedPointSetDimension);

  for (FixedPointIterator it = m_FixedPointSet->GetPoints()->Begin(); it != m_FixedPointSet->GetPoints()->End(); ++it) {
	  typename FixedPointSetType::PointType point = it.Value();

    for (size_t dim = 0; dim < FixedPointSetDimension; ++dim) {
      m_FixedPointMatrix(it.Index(), dim) = point[dim];
    }
  }

  std::cout << "points in scene " << m_FixedPointMatrix.cols() << " " << m_FixedPointMatrix.rows() << std::endl;

  // If the PointSet is provided by a source, update the source.
  if ( m_MovingPointSet->GetSource() )
    {
	  m_MovingPointSet->GetSource()->Update();
    }

  m_MovingPointMatrix.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);
  m_TransformedPointMatrix.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    typename MovingPointSetType::PointType point = it.Value();

    for (size_t dim = 0; dim < MovingPointSetDimension; ++dim) {
      m_MovingPointMatrix(it.Index(), dim) = point[dim];
    }
  }

  m_Gradient.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);
}

/** PrintSelf */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving PointSet: " << m_MovingPointSet.GetPointer()  << std::endl;
  os << indent << "Fixed  PointSet: " << m_FixedPointSet.GetPointer()   << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
}
} // end namespace itk

#endif
