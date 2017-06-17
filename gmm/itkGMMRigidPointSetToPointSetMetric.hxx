#ifndef itkGMMRigidPointSetToPointSetMetric_hxx
#define itkGMMRigidPointSetToPointSetMetric_hxx

#include "itkGMMRigidPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMRigidPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMRigidPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize(void)
throw (ExceptionObject)
{
  Superclass::Initialize();
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}

/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void GMMRigidPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  m_Transform->SetParameters(parameters);

  if (derivative.size() != m_NumberOfParameters) {
    derivative.set_size(m_NumberOfParameters);
  }

  // compute transformed point set
  m_TransformedPointSet = MovingPointSetType::New();
  for (MovingPointIterator iter = m_MovingPointSet->GetPoints()->Begin(); iter != m_MovingPointSet->GetPoints()->End(); ++iter) {
    m_TransformedPointSet->SetPoint(iter.Index(), m_Transform->TransformPoint(iter.Value()));
  }

  GradientType gradient;
  gradient.Fill(0);

  double scale = (m_FixedPointSetScale*m_FixedPointSetScale + m_MovingPointSetScale*m_MovingPointSetScale) / 2;

  for (MovingPointIterator movingIter = m_TransformedPointSet->GetPoints()->Begin(); movingIter != m_TransformedPointSet->GetPoints()->End(); ++movingIter) {
    const typename MovingPointSetType::PointType transformedPoint = movingIter.Value();

    for (FixedPointIterator fixedIter = m_FixedPointSet->GetPoints()->Begin(); fixedIter != m_FixedPointSet->GetPoints()->End(); ++fixedIter) {
      const typename FixedPointSetType::PointType fixedPoint = fixedIter.Value();
      const double distance = transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);
      const double expval = exp(-distance / scale);
      value += expval;

      for (size_t dim = 0; dim < PointDimension; ++dim) {
        gradient[dim] += (-2.0) * expval * (transformedPoint[dim] - fixedPoint[dim]) / scale;
      }
    }

    // compute the derivatives
    m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(movingIter.Value(), m_Jacobian, m_JacobianCache);

    for (size_t par = 0; par < m_NumberOfParameters; par++) {
      for (size_t dim = 0; dim < Self::FixedPointSetDimension; dim++) {
        derivative[par] += m_Jacobian(dim, par) * gradient[dim];
      }
    }
  }

  const double factor = m_TransformedPointSet->GetNumberOfPoints() * m_FixedPointSet->GetNumberOfPoints();

  value /= factor;
  value *= -1.0;

  for (size_t par = 0; par < m_NumberOfParameters; par++) {
    derivative[par] = (-2.0) * derivative[par] / factor;
  }
}
}

#endif
