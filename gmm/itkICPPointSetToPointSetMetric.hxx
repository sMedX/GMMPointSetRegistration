#ifndef itkICPPointSetToPointSetMetric_hxx
#define itkICPPointSetToPointSetMetric_hxx

#include "itkICPPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::ICPPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
ICPPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  typename FixedPointSetType::PointsContainer::ConstPointer fixedPointContainer = this->m_FixedPointSet->GetPoints();

  m_FixedPointsLocator = FixedPointsLocatorType::New();
  m_FixedPointsLocator->SetPoints(const_cast<typename FixedPointSetType::PointsContainer*> (fixedPointContainer.GetPointer()));
  m_FixedPointsLocator->Initialize();
}

/**
 * Get the match Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
typename ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::MeasureType
ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValue(const TransformParametersType & parameters) const
{
  itkExceptionMacro(<< "not implemented");
}
/**
 * Get the Derivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
 * Get both the match Measure and theDerivative Measure
 */
template <typename TFixedPointSet, typename TMovingPointSet>
void ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  this->m_Transform->SetParameters(parameters);

  if (derivative.size() != this->m_NumberOfParameters) {
    derivative.set_size(this->m_NumberOfParameters);
  }

  derivative.Fill(NumericTraits<typename DerivativeType::ValueType>::ZeroValue());
  value = NumericTraits<MeasureType>::ZeroValue();
  GradientType gradient;

  for (MovingPointIterator movingIter = this->m_MovingPointSet->GetPoints()->Begin(); movingIter != this->m_MovingPointSet->GetPoints()->End(); ++movingIter) {
    const typename MovingPointSetType::PointType movingPoint = movingIter.Value();
    const typename MovingPointSetType::PointType transformedPoint = this->m_Transform->TransformPoint(movingPoint);

    // find closest point
    size_t idx = m_FixedPointsLocator->FindClosestPoint(transformedPoint);
    const typename FixedPointSetType::PointType fixedPoint = this->m_FixedPointSet->GetPoint(idx);

    // compute value
    value += transformedPoint.SquaredEuclideanDistanceTo(fixedPoint);

    // compute gradient for the current moving point
    for (size_t dim = 0; dim < PointDimension; ++dim) {
      gradient[dim] = 2.0 * (transformedPoint[dim] - fixedPoint[dim]);
    }

    // compute derivatives
    this->m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(movingPoint, this->m_Jacobian, this->m_JacobianCache);

    for (size_t par = 0; par < this->m_NumberOfParameters; ++par) {
      for (size_t dim = 0; dim < PointDimension; ++dim) {
        derivative[par] += this->m_Jacobian(dim, par) * gradient[dim];
      }
    }
  }

  value /= this->m_MovingPointSet->GetNumberOfPoints();

  for (size_t par = 0; par < this->m_NumberOfParameters; ++par) {
    derivative[par] /= this->m_MovingPointSet->GetNumberOfPoints();
  }
}
}

#endif
