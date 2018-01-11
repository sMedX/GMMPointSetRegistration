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

  m_TransformedMovingPointSet = ITK_NULLPTR;

  m_Jacobian.set_size(MovingPointSetDimension, m_NumberOfParameters);
  m_JacobianCache.set_size(MovingPointSetDimension, MovingPointSetDimension);

  m_NormalizingValueFactor = 1;
  m_NormalizingDerivativeFactor = 1;
  m_Scale = 1;

  m_UseFixedPointSetKdTree = false;
  m_FixedPointsLocator = ITK_NULLPTR;

  m_UseMovingPointSetKdTree = false;
  m_MovingPointsLocator = ITK_NULLPTR;

  m_Radius = 3;
}

/**
* Get the match Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>::MeasureType
GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  this->InitializeForIteration(parameters);

  MeasureType value = NumericTraits<MeasureType>::ZeroValue();

  for (MovingPointIterator it = m_TransformedMovingPointSet->GetPoints()->Begin(); it != m_TransformedMovingPointSet->GetPoints()->End(); ++it) 
  {
    value += GetLocalNeighborhoodValue(it.Value());
  }

  value *= m_NormalizingValueFactor;

  return value;
}

/*
* Get both the match Measure and the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void 
GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>
::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->InitializeForIteration(parameters);

  value = NumericTraits<MeasureType>::ZeroValue();
  MeasureType localValue;

  if (derivative.size() != this->m_NumberOfParameters) 
  {
    derivative.set_size(this->m_NumberOfParameters);
  }

  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());
  LocalDerivativeType localDerivative;

  for (MovingPointIterator it = m_TransformedMovingPointSet->GetPoints()->Begin(); it != m_TransformedMovingPointSet->GetPoints()->End(); ++it) 
  {
    // compute local value and derivatives
    this->GetLocalNeighborhoodValueAndDerivative(it.Value(), localValue, localDerivative);

    value += localValue;

    // compute derivatives
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(m_MovingPointSet->GetPoint(it.Index()), m_Jacobian, m_JacobianCache);

    for (size_t dim = 0; dim < PointDimension; ++dim) 
    {
      for (size_t par = 0; par < m_NumberOfParameters; ++par) 
      {
        derivative[par] += m_Jacobian(dim, par) * localDerivative[dim];
      }
    }
  }

  value *= m_NormalizingValueFactor;

  for (size_t par = 0; par < m_NumberOfParameters; ++par) 
  {
    derivative[par] *= m_NormalizingDerivativeFactor;
  }
}

/**
* Get the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/** Initialize data for current iteration with the input parameters */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::InitializeForIteration(const ParametersType & parameters) const
{
  this->SetTransformParameters(parameters);

  if (!m_TransformedMovingPointSet) 
  {
    m_TransformedMovingPointSet = MovingPointSetType::New();
    m_TransformedMovingPointSet->GetPoints()->resize(m_MovingPointSet->GetNumberOfPoints());
  }

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) 
  {
    m_TransformedMovingPointSet->GetPoints()->SetElement(it.Index(), m_Transform->TransformPoint(it.Value()));
  }
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
  if (!m_Transform) 
    {
    itkExceptionMacro(<< "Transform is not present");
    }

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

  // If the PointSet is provided by a source, update the source.
  if ( m_MovingPointSet->GetSource() )
    {
	  m_MovingPointSet->GetSource()->Update();
    }

  // initialize KdTrees 
  if (m_UseFixedPointSetKdTree && !m_FixedPointsLocator)
    {
    InitializeFixedTree();
    }

  if (m_UseMovingPointSetKdTree && !m_MovingPointsLocator)
    {
    InitializeMovingTree();
  }

  m_NumberOfParameters = m_Transform->GetNumberOfParameters();
  m_NumberOfFixedPoints = m_FixedPointSet->GetNumberOfPoints();
  m_NumberOfMovingPoints = m_MovingPointSet->GetNumberOfPoints();
}

/** Initialize KdTree for FixedPointSet */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::InitializeFixedTree()
{
  m_FixedPointsLocator = FixedPointsLocatorType::New();
  m_FixedPointsLocator->SetPoints(const_cast<FixedPointsContainer*>(m_FixedPointSet->GetPoints()));
  m_FixedPointsLocator->Initialize();
}

/** Initialize KdTree for MovingPointSet */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::InitializeMovingTree()
{
  m_MovingPointsLocator = MovingPointsLocatorType::New();
  m_MovingPointsLocator->SetPoints(const_cast<MovingPointsContainer*> (m_MovingPointSet->GetPoints()));
  m_MovingPointsLocator->Initialize();
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
  os << indent << "Transform:       " << m_Transform.GetPointer()    << std::endl;
}
} // end namespace itk

#endif
