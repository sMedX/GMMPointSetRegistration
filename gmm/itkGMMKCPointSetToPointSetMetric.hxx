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
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();

  const double factor = (double) this->m_MovingPointSet->GetNumberOfPoints() / (double) this->m_FixedPointSet->GetNumberOfPoints();

  this->m_NormalizingValueFactor = -factor / (this->m_MovingPointSet->GetNumberOfPoints() * this->m_FixedPointSet->GetNumberOfPoints());

  this->m_NormalizingDerivativeFactor = -4.0 * factor * this->m_NormalizingValueFactor;
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointType & point) const
{
  const double scale = this->m_Scale*this->m_Scale;

  // compute value for the first sum
  double value1 = 0;

  for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale);
    value1 += expval;
  }

  // compute value for the second sum
  double value2 = 0;

  for (MovingPointIterator it = this->m_TransformedPointSet->GetPoints()->Begin(); it != this->m_TransformedPointSet->GetPoints()->End(); ++it) {
    const double distance = point.SquaredEuclideanDistanceTo(it.Value());
    const double expval = std::exp(-distance / scale);
    value2 += expval;
  }

  // compute local value
  const double ratio = value1 / value2;
  const double value = value1 * ratio;

  return value;
}

template<typename TFixedPointSet, typename TMovingPointSet>
void
GMMKCPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointType & point, MeasureType & value, LocalDerivativeType & derivative) const
{
  const double scale = this->m_Scale * this->m_Scale;

  // compute gradient for the first sum
  double value1 = 0;
  LocalDerivativeType derivative1;
  derivative1.Fill(0);

  for (FixedPointIterator it = this->m_FixedPointSet->GetPoints()->Begin(); it != this->m_FixedPointSet->GetPoints()->End(); ++it) {
    const typename FixedPointSetType::PointType fixedPoint = it.Value();
    const double distance = point.SquaredEuclideanDistanceTo(fixedPoint);
    const double expval = std::exp(-distance / scale);
    value1 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative1[dim] += expval * (point[dim] - fixedPoint[dim]) / scale;
    }
  }

  // compute gradient for the second part
  double value2 = 0;
  LocalDerivativeType derivative2;
  derivative2.Fill(0);

  for (MovingPointIterator it = this->m_TransformedPointSet->GetPoints()->Begin(); it != this->m_TransformedPointSet->GetPoints()->End(); ++it) {
    const typename MovingPointSetType::PointType transformedPoint = it.Value();
    const double distance = point.SquaredEuclideanDistanceTo(transformedPoint);
    const double expval = std::exp(-distance / scale);
    value2 += expval;

    for (size_t dim = 0; dim < this->PointDimension; ++dim) {
      derivative2[dim] += expval * (point[dim] - transformedPoint[dim]) / scale;
    }
  }

  // compute local value
  const double ratio = value1 / value2;
  value = value1 * ratio;

  // compute local derivatives
  for (size_t dim = 0; dim < this->PointDimension; ++dim) {
    derivative[dim] = (derivative1[dim] - derivative2[dim] * ratio) * ratio;
  }
}
}

#endif
