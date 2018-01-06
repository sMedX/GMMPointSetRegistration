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

  this->m_NormalizingValueFactor = 1.0 / this->m_MovingPointSet->GetNumberOfPoints();
}

template<typename TFixedPointSet, typename TMovingPointSet>
typename ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::MeasureType
ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValue(const MovingPointType & point) const
{
  // find closest point
  size_t idx = m_FixedPointsLocator->FindClosestPoint(point);
  const typename FixedPointSetType::PointType fixedPoint = this->m_FixedPointSet->GetPoint(idx);

  return point.SquaredEuclideanDistanceTo(fixedPoint);
}

template<typename TFixedPointSet, typename TMovingPointSet>
void
ICPPointSetToPointSetMetric<TFixedPointSet, TMovingPointSet>
::GetLocalNeighborhoodValueAndDerivative(const MovingPointType & point, MeasureType & value, LocalDerivativeType & derivative) const
{
  // find closest point
  size_t idx = m_FixedPointsLocator->FindClosestPoint(point);
  const typename FixedPointSetType::PointType fixedPoint = this->m_FixedPointSet->GetPoint(idx);

  // compute value
  value = point.SquaredEuclideanDistanceTo(fixedPoint);

  // compute gradient for the current moving point
  for (size_t dim = 0; dim < this->PointDimension; ++dim) {
    derivative[dim] = 2.0 * (point[dim] - fixedPoint[dim]);
  }
}
}

#endif
