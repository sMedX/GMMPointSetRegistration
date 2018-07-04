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

  m_Scale = 1;

  m_UseFixedPointSetKdTree = false;
  m_FixedPointsLocator = ITK_NULLPTR;

  m_UseMovingPointSetKdTree = false;
  m_MovingPointsLocator = ITK_NULLPTR;

  m_SearchRadius = 3;
  m_IsInitialized = false;
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

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) 
  {
    value += GetLocalNeighborhoodValue(it);
  }

  double valueFactor = this->GetNormalizingValueFactor();

  value *= valueFactor;

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

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) 
  {
    // compute local value and derivatives
    if (this->GetLocalNeighborhoodValueAndDerivative(it, localValue, localDerivative))
    {
      value += localValue;

      // compute derivatives
      this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(it.Value(), m_Jacobian, m_JacobianCache);

      for (size_t dim = 0; dim < PointDimension; ++dim) 
      {
        for (size_t par = 0; par < m_NumberOfParameters; ++par) 
        {
          derivative[par] += m_Jacobian(dim, par) * localDerivative[dim];
        }
      }
    }
  }

  double valueFactor = this->GetNormalizingValueFactor();
  double derivativeFactor = this->GetNormalizingDerivativeFactor();

  value *= valueFactor;

  for (size_t par = 0; par < m_NumberOfParameters; ++par) 
  {
    derivative[par] *= derivativeFactor;
  }
}

/*
* Get both the match Measure and the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void
GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>
::GetDerivatives(const TransformParametersType & parameters, DerivativeType & derivative1, DerivativeType & derivative2) const
{
  if (derivative1.size() != this->m_NumberOfParameters) 
  {
    derivative1.set_size(this->m_NumberOfParameters);
  }

  derivative1.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  if (derivative2.size() != this->m_NumberOfParameters) 
  {
    derivative2.set_size(this->m_NumberOfParameters);
  }

  derivative2.Fill(NumericTraits<DerivativeValueType>::ZeroValue());

  LocalDerivativeType localDerivative1;
  LocalDerivativeType localDerivative2;

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {

    // compute local derivatives
    if (this->GetLocalNeighborhoodDerivatives(it, localDerivative1, localDerivative2)) {

      // compute derivatives
      this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(it.Value(), m_Jacobian, m_JacobianCache);

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        for (size_t par = 0; par < m_NumberOfParameters; ++par) {
          derivative1[par] += m_Jacobian(dim, par) * localDerivative1[dim];
          derivative2[par] += m_Jacobian(dim, par) * localDerivative2[dim];
        }
      }
    }
  }

  double derivativeFactor = this->GetNormalizingDerivativeFactor();
  double derivativeScaleFactor = this->GetNormalizingDerivativeScaleFactor();

  for (size_t par = 0; par < m_NumberOfParameters; ++par) 
  {
    derivative1[par] *= derivativeFactor;
    derivative2[par] *= derivativeScaleFactor;
  }
}

/**
* Get the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  this->InitializeForIteration(parameters);

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }

  derivative.Fill(NumericTraits<DerivativeValueType>::ZeroValue());
  LocalDerivativeType localDerivative;

  for (MovingPointIterator it = m_MovingPointSet->GetPoints()->Begin(); it != m_MovingPointSet->GetPoints()->End(); ++it) {
    // compute local value and derivatives
    if (this->GetLocalNeighborhoodDerivative(it, localDerivative)) {

      // compute derivatives
      this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(it.Value(), m_Jacobian, m_JacobianCache);

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        for (size_t par = 0; par < m_NumberOfParameters; ++par) {
          derivative[par] += m_Jacobian(dim, par) * localDerivative[dim];
        }
      }
    }
  }

  double derivativeFactor = this->GetNormalizingDerivativeFactor();

  for (size_t par = 0; par < m_NumberOfParameters; ++par) {
    derivative[par] *= derivativeFactor;
  }
}

/** Set scale */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::SetScale(const double & sigma)
{
  this->m_Scale = sigma; 
  this->m_Variance = sigma * sigma;
}

/** Set scale */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::SetScale(const ParametersType & parameters)
{
  this->SetScale(parameters[0]);
}

/** Initialize data for current iteration with the input parameters */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::InitializeForIteration(const ParametersType & parameters) const
{
  m_Transform->SetParameters(parameters);
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

  m_NumberOfParameters = this->GetNumberOfParameters();

  m_IsInitialized = true;
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

/** Initialize transformed point set*/
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::InitializeTransformedMovingPointSet()
{
  m_TransformedMovingPointSet = FixedPointSetType::New();
  m_TransformedMovingPointSet->GetPoints()->resize(m_MovingPointSet->GetNumberOfPoints());
}

/** PrintSelf */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << std::endl;

  os << "Transform " << m_Transform.GetPointer() << std::endl;
  os << indent << m_Transform->GetNameOfClass() << std::endl;
  os << std::endl;

  os << "Fixed point set  " << m_FixedPointSet.GetPointer() << std::endl;
  os << indent << "Use KdTree  " << m_UseFixedPointSetKdTree << std::endl;
  os << std::endl;

  os << "Moving point set " << m_MovingPointSet.GetPointer() << std::endl;
  os << indent << "Use KdTree  " << m_UseMovingPointSetKdTree << std::endl;
  os << std::endl;
}
} // end namespace itk

#endif
