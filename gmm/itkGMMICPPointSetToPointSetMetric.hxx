#ifndef itkGMMICPPointSetToPointSetMetric_hxx
#define itkGMMICPPointSetToPointSetMetric_hxx

#include "itkGMMICPPointSetToPointSetMetric.h"

namespace itk
{
/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingPointSet>
GMMICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>::GMMICPPointSetToPointSetMetric()
{
}

/** Initialize the metric */
template< typename TFixedPointSet, typename TMovingPointSet >
void
GMMICPPointSetToPointSetMetric< TFixedPointSet, TMovingPointSet >
::Initialize() throw (ExceptionObject)
{
  Superclass::Initialize();
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename GMMICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
GMMICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointIterator & it) const
{
  // transform point
  const FixedPointType point = this->m_Transform->TransformPoint(it.Value());

  // find closest point
  size_t idx = this->m_FixedPointsLocator->FindClosestPoint(point);
  const FixedPointType fixedPoint = this->GetFixedPoint(idx);

  return point.SquaredEuclideanDistanceTo(fixedPoint);
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMMICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointIterator & it, MeasureType & value, LocalDerivativeType & derivative) const
{
  // transform point
  const FixedPointType point = this->m_Transform->TransformPoint(it.Value());

  // find closest point
  size_t idx = this->m_FixedPointsLocator->FindClosestPoint(point);
  const FixedPointType fixedPoint = this->GetFixedPoint(idx);

  // compute value
  value = point.SquaredEuclideanDistanceTo(fixedPoint);

  // compute gradient for the current moving point
  for (size_t dim = 0; dim < this->PointDimension; ++dim) {
    derivative[dim] = point[dim] - fixedPoint[dim];
  }

  return true;
}

template<typename TFixedPointSet, typename TMovingPointSet>
bool
GMMICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodDerivative(const MovingPointIterator & it, LocalDerivativeType & derivative) const
{
  // transform point
  const FixedPointType point = this->m_Transform->TransformPoint(it.Value());

  // find closest point
  size_t idx = this->m_FixedPointsLocator->FindClosestPoint(point);
  const FixedPointType fixedPoint = this->GetFixedPoint(idx);

  // compute gradient for the current moving point
  for (size_t dim = 0; dim < this->PointDimension; ++dim) {
    derivative[dim] = point[dim] - fixedPoint[dim];
  }

  return true;
}

}

#endif
