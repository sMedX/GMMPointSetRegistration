#ifndef itkGMMKCPointSetToPointSetMetric_hxx
#define itkGMMKCPointSetToPointSetMetric_hxx

#include "itkGMMKCPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMKCPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMKCPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void)
throw (ExceptionObject)
{
  Superclass::Initialize();

  m_Gradient1.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);
  m_Gradient2.set_size(m_MovingPointSet->GetNumberOfPoints(), MovingPointSetDimension);
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}
/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  m_Transform->SetParameters(parameters);

  if (derivative.size() != m_NumberOfParameters) {
    derivative.set_size(m_NumberOfParameters);
  }

  m_TransformedPointSet = MovingPointSetType::New();
  for (MovingPointIterator iter = m_MovingPointSet->GetPoints()->Begin(); iter != m_MovingPointSet->GetPoints()->End(); ++iter) {
    m_TransformedPointSet->SetPoint(iter.Index(), m_Transform->TransformPoint(iter.Value()));
  }

  double value1 = 0;
  double value2 = 0;

  GradientType gradient1;
  GradientType gradient2;

  LocalDerivativeType derivative1(m_NumberOfParameters);
  derivative1.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  LocalDerivativeType derivative2(m_NumberOfParameters);
  derivative2.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (MovingPointIterator movingIter1 = m_TransformedPointSet->GetPoints()->Begin(); movingIter1 != m_TransformedPointSet->GetPoints()->End(); ++movingIter1) {
    const typename MovingPointSetType::PointType transformedPoint1 = movingIter1.Value();

    //------------------------------------
    // compute gradient for the first part
    gradient1.Fill(0);

    for (FixedPointIterator fixedIter = m_FixedPointSet->GetPoints()->Begin(); fixedIter != m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
      const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();

      double distance = 0;
      for (size_t dim = 0; dim < PointDimension; ++dim) {
        distance += pow(transformedPoint1[dim] / m_MovingPointSetScale - fixedPoint[dim] / m_FixedPointSetScale, 2);
      }
      
      value1 += exp(-distance);

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        gradient1[dim] += (-2.0) * exp(-distance) * (transformedPoint1[dim] / m_MovingPointSetScale - fixedPoint[dim] / m_FixedPointSetScale) / m_MovingPointSetScale;
      }
    }

    //-------------------------------------
    // compute gradient for the second part
    gradient2.Fill(0);

    for (MovingPointIterator movingIter2 = m_TransformedPointSet->GetPoints()->Begin(); movingIter2 != m_TransformedPointSet->GetPoints()->End(); ++movingIter2) {
      const typename MovingPointSetType::PointType transformedPoint2 = movingIter2.Value();

      double distance = 0;
      for (size_t dim = 0; dim < PointDimension; ++dim) {
        distance += pow(transformedPoint1[dim] / m_MovingPointSetScale - transformedPoint2[dim] / m_MovingPointSetScale, 2);
      }

      value2 += exp(-distance);

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        gradient2[dim] += (-2.0) * exp(-distance) * (transformedPoint1[dim] / m_MovingPointSetScale - transformedPoint2[dim] / m_MovingPointSetScale) / m_MovingPointSetScale;
      }
    }

    //-------------------------------------
    // compute derivatives
    m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(m_MovingPointSet->GetPoint(movingIter1.Index()), m_Jacobian, m_JacobianCache);

    for (size_t par = 0; par < m_NumberOfParameters; par++) {
      for (size_t dim = 0; dim < PointDimension; dim++) {
        derivative1[par] += m_Jacobian(dim, par) * gradient1[dim];
        derivative2[par] += m_Jacobian(dim, par) * gradient2[dim];
      }
    }
  }

  const double factor = m_MovingPointSet->GetNumberOfPoints() * m_FixedPointSet->GetNumberOfPoints();
  const double ratio = value1 / value2;
  
  value = -value1 * ratio / factor;

  for (size_t par = 0; par < m_NumberOfParameters; par++) {
    derivative[par] = (-2.0) * (derivative1[par] - derivative2[par] * ratio) * ratio / factor;
  }
}
}

#endif
