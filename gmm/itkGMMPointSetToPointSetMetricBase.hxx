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

  m_NormalizingValueFactor = 1;

  m_UseKdTree = true;
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

  for (MovingPointIterator it = this->m_TransformedPointSet->GetPoints()->Begin(); it != this->m_TransformedPointSet->GetPoints()->End(); ++it) {
    value += GetLocalNeighborhoodValue(it.Value());
  }

  value *= m_NormalizingValueFactor;

  return value;
}

/*
* Get both the match Measure and the Derivative Measure
*/
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMPointSetToPointSetMetricBase<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->InitializeForIteration(parameters);

  double scale = this->m_Scale * this->m_Scale;

  value = NumericTraits<MeasureType>::ZeroValue();

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }
  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  MeasureType localValue;

  LocalDerivativeType localDerivative;

  for (MovingPointIterator it = this->m_TransformedPointSet->GetPoints()->Begin(); it != this->m_TransformedPointSet->GetPoints()->End(); ++it) {

    // compute local value and derivatives
    this->GetLocalNeighborhoodValueAndDerivative(it.Value(), localValue, localDerivative);

    value += localValue;

    // compute derivatives
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(this->m_MovingPointSet->GetPoint(it.Index()), this->m_Jacobian, this->m_JacobianCache);

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      for (size_t par = 0; par < this->m_NumberOfParameters; ++par) {
        derivative[par] += this->m_Jacobian(dim, par) * localDerivative[dim];
      }
    }
  }

  value *= m_NormalizingValueFactor;

  for (size_t par = 0; par < this->m_NumberOfParameters; par++) {
    derivative[par] *= m_NormalizingValueFactor;
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

/** Set the parameters that define a unique transform */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::InitializeForIteration(const ParametersType & parameters) const
{
  this->SetTransformParameters(parameters);

  this->m_TransformedPointSet = MovingPointSetType::New();

  for (MovingPointIterator iter = this->m_MovingPointSet->GetPoints()->Begin(); iter != this->m_MovingPointSet->GetPoints()->End(); ++iter) {
    this->m_TransformedPointSet->SetPoint(iter.Index(), this->m_Transform->TransformPoint(iter.Value()));
  }
}

/** Set the parameters that define a unique transform */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMPointSetToPointSetMetricBase< TFixedPointSet, TMovingPointSet >
::SetTransformParameters(const ParametersType & parameters) const
{
  if ( !m_Transform ) {
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

  // If the PointSet is provided by a source, update the source.
  if ( m_MovingPointSet->GetSource() )
    {
	  m_MovingPointSet->GetSource()->Update();
    }

  m_NumberOfFixedPoints = m_FixedPointSet->GetNumberOfPoints();
  m_NumberOfMovingPoints = m_MovingPointSet->GetNumberOfPoints();
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
